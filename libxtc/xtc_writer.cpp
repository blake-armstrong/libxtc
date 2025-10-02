/*
 * Copyright (c) 2025
 * All rights reserved.
 *
 * C++ XTCWriter class - implementation file
 */

#include "xtc_writer.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

// Endian swap utilities
template <typename T> T XTCWriter::byteswap(T value) {
  union {
    T val;
    uint8_t bytes[sizeof(T)];
  } src, dst;

  src.val = value;
  for (size_t i = 0; i < sizeof(T); i++) {
    dst.bytes[i] = src.bytes[sizeof(T) - 1 - i];
  }
  return dst.val;
}

// XDR packing (big-endian format)
void XTCWriter::pack_int(std::vector<uint8_t> &buf, int32_t value) {
  int32_t swapped = byteswap(value);
  const uint8_t *ptr = reinterpret_cast<const uint8_t *>(&swapped);
  buf.insert(buf.end(), ptr, ptr + sizeof(int32_t));
}

void XTCWriter::pack_uint(std::vector<uint8_t> &buf, uint32_t value) {
  uint32_t swapped = byteswap(value);
  const uint8_t *ptr = reinterpret_cast<const uint8_t *>(&swapped);
  buf.insert(buf.end(), ptr, ptr + sizeof(uint32_t));
}

void XTCWriter::pack_float(std::vector<uint8_t> &buf, float value) {
  uint32_t tmp;
  std::memcpy(&tmp, &value, sizeof(float));
  tmp = byteswap(tmp);
  const uint8_t *ptr = reinterpret_cast<const uint8_t *>(&tmp);
  buf.insert(buf.end(), ptr, ptr + sizeof(float));
}

XTCWriter::XTCWriter(const std::string &fname, float precision)
    : filename(fname), precision(precision) {

  if (precision <= 0) {
    throw std::invalid_argument("Precision must be positive");
  }

  // Open file for writing
  fxtc.open(fname, std::ios::binary | std::ios::out | std::ios::trunc);
  if (!fxtc.is_open()) {
    throw std::runtime_error("Cannot open file for writing: " + fname);
  }
}

XTCWriter::~XTCWriter() { close(); }

void XTCWriter::close() {
  if (fxtc.is_open()) {
    fxtc.close();
  }
}

namespace {
// Helper functions for compression

int sizeofint(unsigned int size) {
  unsigned int num = 1;
  unsigned int num_of_bits = 0;

  while (size >= num && num_of_bits < 32) {
    num_of_bits++;
    num <<= 1;
  }
  return num_of_bits;
}

int sizeofints(int num_of_ints, unsigned int sizes[]) {
  int i, num;
  unsigned int num_of_bytes, num_of_bits, bytes[32], bytecnt, tmp;
  num_of_bytes = 1;
  bytes[0] = 1;
  num_of_bits = 0;
  for (i = 0; i < num_of_ints; i++) {
    tmp = 0;
    for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++) {
      tmp = bytes[bytecnt] * sizes[i] + tmp;
      bytes[bytecnt] = tmp & 0xff;
      tmp >>= 8;
    }
    while (tmp != 0) {
      bytes[bytecnt++] = tmp & 0xff;
      tmp >>= 8;
    }
    num_of_bytes = bytecnt;
  }
  num = 1;
  num_of_bytes--;
  while (static_cast<int>(bytes[num_of_bytes]) >= num) {
    num_of_bits++;
    num *= 2;
  }
  return num_of_bits + num_of_bytes * 8;
}

static const int magicints[] = {
    0,       0,       0,       0,       0,       0,       0,        0,
    0,       8,       10,      12,      16,      20,      25,       32,
    40,      50,      64,      80,      101,     128,     161,      203,
    256,     322,     406,     512,     645,     812,     1024,     1290,
    1625,    2048,    2580,    3250,    4096,    5060,    6501,     8192,
    10321,   13003,   16384,   20642,   26007,   32768,   41285,    52015,
    65536,   82570,   104031,  131072,  165140,  208063,  262144,   330280,
    416127,  524287,  660561,  832255,  1048576, 1321122, 1664510,  2097152,
    2642245, 3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085,
    16777216};

const int FIRSTIDX = 9;
const int LASTIDX = sizeof(magicints) / sizeof(*magicints);

class bit_writer {
public:
  bit_writer() : buf(0), bitsused(0) {}

  std::vector<uint8_t> get_buffer() const {
    std::vector<uint8_t> result = written;
    if (bitsused > 0) {
      // Flush remaining bits
      int bytes_left = (bitsused + 7) / 8;
      uint64_t tmp = buf;
      for (int i = 0; i < bytes_left; i++) {
        result.push_back((tmp >> 56) & 0xFF);
        tmp <<= 8;
      }
    }
    return result;
  }

  void write_bits(uint64_t value, int nbits) {
    buf |= (value << (64 - bitsused - nbits));
    bitsused += nbits;

    while (bitsused >= 8) {
      written.push_back((buf >> 56) & 0xFF);
      buf <<= 8;
      bitsused -= 8;
    }
  }

  void write_int(int value, int nbits) {
    uint64_t mask = (1ULL << nbits) - 1;
    write_bits(value & mask, nbits);
  }

private:
  uint64_t buf;
  int bitsused;
  std::vector<uint8_t> written;
};

void encodeints(bit_writer &bw, int num_of_ints, int num_of_bits,
                unsigned int sizes[], int nums[]) {
  int bytes[32];
  int i, j, num_of_bytes;

  unsigned int tmp = static_cast<unsigned int>(nums[0]);
  bytes[0] = tmp & 0xff;
  bytes[1] = (tmp >> 8) & 0xff;
  bytes[2] = (tmp >> 16) & 0xff;
  bytes[3] = (tmp >> 24) & 0xff;

  num_of_bytes = 4;
  for (i = 1; i < num_of_ints; i++) {
    if (nums[i] >= static_cast<int>(sizes[i])) {
      throw std::runtime_error("Compression error: value out of bounds");
    }

    unsigned int carry = nums[i];
    for (j = 0; j < num_of_bytes; j++) {
      tmp = bytes[j] * sizes[i] + carry;
      bytes[j] = tmp & 0xff;
      carry = tmp >> 8;
    }
    while (carry != 0) {
      bytes[num_of_bytes++] = carry & 0xff;
      carry >>= 8;
    }
  }

  // Write bytes to bit stream
  int bits_left = num_of_bits;
  for (i = 0; i < num_of_bytes && bits_left > 0; i++) {
    int bits_to_write = std::min(8, bits_left);
    bw.write_bits(bytes[i], bits_to_write);
    bits_left -= bits_to_write;
  }
}

} // namespace

void XTCWriter::calculate_bounds(uint32_t natoms, const float *coords,
                                 int32_t minint[3], int32_t maxint[3],
                                 std::vector<int32_t> &int_coords) {
  int_coords.resize(natoms * 3);

  // Initialize bounds
  for (int d = 0; d < 3; d++) {
    minint[d] = std::numeric_limits<int32_t>::max();
    maxint[d] = std::numeric_limits<int32_t>::min();
  }

  // Convert to integer coordinates and find bounds
  for (uint32_t i = 0; i < natoms; i++) {
    for (int d = 0; d < 3; d++) {
      int32_t ival =
          static_cast<int32_t>(std::round(coords[i * 3 + d] * precision));
      int_coords[i * 3 + d] = ival;
      minint[d] = std::min(minint[d], ival);
      maxint[d] = std::max(maxint[d], ival);
    }
  }
}

bool XTCWriter::compress_coords(uint32_t natoms, const float *coords,
                                std::vector<uint8_t> &output, int32_t &nbytes) {
  int32_t minint[3], maxint[3];
  calculate_bounds(natoms, coords, minint, maxint, int_coords);

  uint32_t sizeint[3];
  for (int i = 0; i < 3; i++) {
    sizeint[i] = maxint[i] - minint[i] + 1;
  }

  // Check if we should use large integer mode
  bool large = (sizeint[0] | sizeint[1] | sizeint[2]) > 0xffffff;
  unsigned int bitsize;
  unsigned int bitsizeint[3];

  if (large) {
    bitsizeint[0] = sizeofint(sizeint[0]);
    bitsizeint[1] = sizeofint(sizeint[1]);
    bitsizeint[2] = sizeofint(sizeint[2]);
    bitsize = 0;
  } else {
    bitsize = sizeofints(3, sizeint);
  }

  bit_writer bw;

  // Determine smallidx
  int smallidx = FIRSTIDX;
  while (smallidx < LASTIDX &&
         magicints[smallidx] < static_cast<int>(sizeint[0]) &&
         magicints[smallidx] < static_cast<int>(sizeint[1]) &&
         magicints[smallidx] < static_cast<int>(sizeint[2])) {
    smallidx++;
  }

  int thiscoord[3], prevcoord[3] = {0, 0, 0};

  for (uint32_t i = 0; i < natoms; i++) {
    for (int d = 0; d < 3; d++) {
      thiscoord[d] = int_coords[i * 3 + d] - minint[d];
    }

    // Write main coordinate
    if (large) {
      for (int d = 0; d < 3; d++) {
        bw.write_int(thiscoord[d], bitsizeint[d]);
      }
    } else {
      encodeints(bw, 3, bitsize, sizeint, thiscoord);
    }

    prevcoord[0] = thiscoord[0];
    prevcoord[1] = thiscoord[1];
    prevcoord[2] = thiscoord[2];

    // Simple implementation: no run-length encoding for now
    // Just write flag = 0 (no small coordinates follow)
    bw.write_bits(0, 1);

    // Update for next iteration (keeping smallidx constant for simplicity)
  }

  output = bw.get_buffer();
  nbytes = output.size();
  return true;
}

bool XTCWriter::write_frame(uint32_t natoms, uint32_t step, float time,
                            const std::array<std::array<float, 3>, 3> &box,
                            const float *coords, bool compress) {
  if (!fxtc.is_open()) {
    return false;
  }

  std::vector<uint8_t> header;

  // Write magic number
  pack_int(header, mxtc);

  // Write natoms
  pack_uint(header, natoms);

  // Write step
  pack_uint(header, step);

  // Write time
  pack_float(header, time);

  // Write box matrix (9 floats)
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      pack_float(header, box[i][j]);
    }
  }

  std::vector<uint8_t> coord_data;
  int32_t nbytes = 0;

  if (compress && natoms > 0) {
    // Prepare compression header
    std::vector<uint8_t> comp_header;
    pack_float(comp_header, precision);

    int32_t minint[3], maxint[3];
    calculate_bounds(natoms, coords, minint, maxint, int_coords);

    for (int i = 0; i < 3; i++) {
      pack_int(comp_header, minint[i]);
    }
    for (int i = 0; i < 3; i++) {
      pack_int(comp_header, maxint[i]);
    }

    uint32_t sizeint[3];
    for (int i = 0; i < 3; i++) {
      sizeint[i] = maxint[i] - minint[i] + 1;
    }

    int smallidx = FIRSTIDX;
    while (smallidx < LASTIDX &&
           magicints[smallidx] < static_cast<int>(sizeint[0]) &&
           magicints[smallidx] < static_cast<int>(sizeint[1]) &&
           magicints[smallidx] < static_cast<int>(sizeint[2])) {
      smallidx++;
    }

    pack_int(comp_header, smallidx);

    // Compress coordinates
    std::vector<uint8_t> compressed;
    compress_coords(natoms, coords, compressed, nbytes);

    pack_int(comp_header, nbytes);

    // Combine headers and data
    uint32_t lsize = comp_header.size() + align4(nbytes);
    pack_uint(header, lsize);

    // Write to file
    fxtc.write(reinterpret_cast<const char *>(header.data()), header.size());
    fxtc.write(reinterpret_cast<const char *>(comp_header.data()),
               comp_header.size());
    fxtc.write(reinterpret_cast<const char *>(compressed.data()), nbytes);

    // Add padding to align to 4 bytes
    int padding = align4(nbytes) - nbytes;
    for (int i = 0; i < padding; i++) {
      uint8_t zero = 0;
      fxtc.write(reinterpret_cast<const char *>(&zero), 1);
    }

  } else {
    // Uncompressed
    uint32_t lsize = natoms * 3 * sizeof(float);
    pack_uint(header, lsize);

    fxtc.write(reinterpret_cast<const char *>(header.data()), header.size());

    // Write coordinates in big-endian format
    for (uint32_t i = 0; i < natoms * 3; i++) {
      uint32_t tmp;
      std::memcpy(&tmp, &coords[i], sizeof(float));
      tmp = byteswap(tmp);
      fxtc.write(reinterpret_cast<const char *>(&tmp), sizeof(float));
    }
  }

  return fxtc.good();
}

bool XTCWriter::write_frame(uint32_t step, float time,
                            const std::array<std::array<float, 3>, 3> &box,
                            const std::vector<float> &coords, bool compress) {
  uint32_t natoms = coords.size() / 3;
  return write_frame(natoms, step, time, box, coords.data(), compress);
}

bool XTCWriter::write_frame(uint32_t step, float time,
                            const std::array<std::array<float, 3>, 3> &box,
                            const std::vector<std::array<float, 3>> &coords,
                            bool compress) {
  std::vector<float> flat_coords(coords.size() * 3);
  for (size_t i = 0; i < coords.size(); i++) {
    flat_coords[i * 3] = coords[i][0];
    flat_coords[i * 3 + 1] = coords[i][1];
    flat_coords[i * 3 + 2] = coords[i][2];
  }
  return write_frame(coords.size(), step, time, box, flat_coords.data(),
                     compress);
}
