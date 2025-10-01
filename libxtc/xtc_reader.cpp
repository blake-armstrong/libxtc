/*
 * Copyright (c) 2020, Nikolay A. Krylov
 * All rights reserved.
 *
 * C++ XTCReader class - implementation file
 */

#include "xtc_reader.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

// Endian swap utilities
template <typename T> T XTCReader::byteswap(T value) {
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

// XDR unpacking (big-endian format)
int32_t XTCReader::unpack_int(const uint8_t *&ptr) {
  int32_t value;
  std::memcpy(&value, ptr, sizeof(int32_t));
  ptr += sizeof(int32_t);
  return byteswap(value);
}

uint32_t XTCReader::unpack_uint(const uint8_t *&ptr) {
  uint32_t value;
  std::memcpy(&value, ptr, sizeof(uint32_t));
  ptr += sizeof(uint32_t);
  return byteswap(value);
}

float XTCReader::unpack_float(const uint8_t *&ptr) {
  uint32_t tmp;
  std::memcpy(&tmp, ptr, sizeof(uint32_t));
  tmp = byteswap(tmp);
  float value;
  std::memcpy(&value, &tmp, sizeof(float));
  ptr += sizeof(float);
  return value;
}

XTCReader::XTCReader(const std::string &fname, float beg, float end, float dt,
                     int nt)
    : t(0.0f), filename(fname), nacur(0), lsize(0), datalen(0), beg(beg),
      end(end), dt(dt), iframe(0), t0(0.0f), dtrrj(0.0f), teps(0.0f),
      bigdt(false) {

  // Initialize frame_data structure
  std::memset(&frd, 0, sizeof(frame_data));
  frd.nt = nt;

  // Open file
  fxtc.open(fname, std::ios::binary | std::ios::in);
  if (!fxtc.is_open()) {
    throw std::runtime_error("Cannot open file: " + fname);
  }

  // Get file size
  fxtc.seekg(0, std::ios::end);
  fsize = fxtc.tellg();
  fxtc.seekg(0, std::ios::beg);

  // Read first frame to determine trajectory parameters
  if (!next_frame(true, true)) {
    throw std::runtime_error("Failed to read first frame");
  }
  this->beg = this->t;
  t0 = this->t;

  // Read second frame to determine dt
  if (next_frame(true, true)) {
    float t1 = this->t;
    dtrrj = t1 - t0;
    teps = 1e-3f * dtrrj;
  }

  // Reset to beginning
  fxtc.clear();
  fxtc.seekg(0, std::ios::beg);

  // Seek to start time if specified
  if (beg > 0) {
    this->beg = beg;
    if (dtrrj > 0) {
      search_t(beg);
    } else {
      skip_to(beg);
    }
  } else {
    this->t = this->beg;
  }

  // Setup dt sampling strategy
  if (dt > 0) {
    if (dt < dtrrj) {
      this->dt = -1;
    } else {
      int scan_cost = 50000;
      int nfr_skip = static_cast<int>(dt / dtrrj);
      int hdr_cost = scan_cost * hdrfull * nfr_skip;
      float datalen_cost = nacur * flt3len * 1.5f / 3.0f;
      bigdt = datalen_cost < hdr_cost;
    }
  }

  iframe = 0;
}

XTCReader::~XTCReader() {
  if (fxtc.is_open()) {
    fxtc.close();
  }
}

bool XTCReader::read_header(bool check) {
  std::vector<uint8_t> mem(hdr1);

  fxtc.read(reinterpret_cast<char *>(mem.data()), hdr1);
  if (fxtc.gcount() != hdr1) {
    return false;
  }

  const uint8_t *ptr = mem.data();

  int32_t magic = unpack_int(ptr);
  if (magic != mxtc) {
    throw std::runtime_error("Bad frame magic!");
  }

  nacur = unpack_uint(ptr);
  uint32_t fr_idx = unpack_uint(ptr);
  (void)fr_idx; // Unused but part of format
  t = unpack_float(ptr);

  // Read box matrix
  for (int i = 0; i < 9; i++) {
    box[i / 3][i % 3] = unpack_float(ptr);
  }

  lsize = unpack_uint(ptr);

  if (check) {
    if (t < 0 || box[0][1] > 0 || box[0][2] > 0) {
      throw std::runtime_error("Bad frame header!");
    }
  }

  return true;
}

void XTCReader::decompress() {
  if (lsize <= 9) {
    // Uncompressed coordinates
    X.resize(nacur * 3);
    fxtc.read(reinterpret_cast<char *>(X.data()), datalen);

    // Convert from big-endian
    for (size_t i = 0; i < X.size(); i++) {
      uint32_t tmp;
      std::memcpy(&tmp, &X[i], sizeof(float));
      tmp = byteswap(tmp);
      std::memcpy(&X[i], &tmp, sizeof(float));
    }
  } else {
    // Compressed coordinates
    int count8 = (datalen + 7) & (~0x7);
    buf.resize(count8);
    fxtc.read(reinterpret_cast<char *>(buf.data()), datalen);

    X.resize(nacur * 3);

    bool ok = unpack_frame(frd, reinterpret_cast<const uint64_t *>(buf.data()),
                           X.data());
    if (!ok) {
      throw std::runtime_error("Bad frame data at time " + std::to_string(t));
    }
  }
}

bool XTCReader::next_frame(bool hdr_only, bool no_check) {
  if (!read_header()) {
    return false;
  }

  // Compute data length
  if (lsize <= 9) {
    datalen = lsize * flt3len;
  } else {
    std::vector<uint8_t> mem2(hdr2);
    fxtc.read(reinterpret_cast<char *>(mem2.data()), hdr2);

    const uint8_t *ptr = mem2.data();
    float prec = unpack_float(ptr);
    frd.inv_p = 1.0f / prec;

    std::array<int32_t, 8> int_arr;
    for (int i = 0; i < 8; i++) {
      int_arr[i] = unpack_int(ptr);
    }

    int nbytes = int_arr[7];
    for (int i = 0; i < 3; i++) {
      frd.minint[i] = int_arr[i];
      frd.maxint[i] = int_arr[i + 3];
    }
    frd.smli = int_arr[6];
    frd.natoms = nacur;

    datalen = align4(nbytes);
  }

  if (hdr_only) {
    fxtc.seekg(datalen, std::ios::cur);
  } else {
    decompress();
  }

  if (no_check) {
    return true;
  }

  float tcur = t;
  if (dt > 0) {
    iframe++;
    float target_t = tnext();

    if (bigdt) {
      search_t(target_t);
    } else {
      skip_to(target_t);
    }
  }

  t = tcur;
  return true;
}

bool XTCReader::eot() {
  if (fxtc.tellg() >= fsize) {
    return true;
  }

  if (end > 0) {
    float target_t = tnext();
    if (target_t - end > teps) {
      return true;
    }
  }

  return false;
}

float XTCReader::tnext() const {
  if (dt > 0) {
    return beg + iframe * dt;
  } else {
    return t + dtrrj;
  }
}

bool XTCReader::search_frame_hdr() {
  std::streampos start_pos = fxtc.tellg();
  std::vector<char> data(static_cast<size_t>(1.5 * (hdrfull + datalen)));
  fxtc.read(data.data(), data.size());
  std::streamsize bytes_read = fxtc.gcount();

  const char magic_bytes[4] = {static_cast<char>((mxtc >> 24) & 0xFF),
                               static_cast<char>((mxtc >> 16) & 0xFF),
                               static_cast<char>((mxtc >> 8) & 0xFF),
                               static_cast<char>(mxtc & 0xFF)};

  size_t offset = 0;
  while (offset < static_cast<size_t>(bytes_read)) {
    auto it = std::search(data.begin() + offset, data.begin() + bytes_read,
                          magic_bytes, magic_bytes + 4);
    if (it == data.begin() + bytes_read) {
      return false;
    }

    size_t magic_pos = std::distance(data.begin(), it);
    std::streampos cur_pos = start_pos + static_cast<std::streamoff>(magic_pos);
    fxtc.clear();
    fxtc.seekg(cur_pos);

    try {
      read_header(true);
      fxtc.seekg(cur_pos);
      return true;
    } catch (...) {
      offset = magic_pos + 4;
    }
  }

  return false;
}

void XTCReader::search_t(float target_t) {
  if (eot()) {
    return;
  }

  std::streamoff off = static_cast<std::streamoff>(
      std::max(0.0f, (hdrfull + datalen) * ((target_t - t0) / dtrrj - 1.5f)));

  if (off >= static_cast<std::streamoff>(fsize)) {
    fxtc.seekg(0, std::ios::end);
    return;
  }

  fxtc.seekg(off);
  if (!search_frame_hdr()) {
    throw std::runtime_error("Failed to locate next frame header!");
  }

  skip_to(target_t);
}

void XTCReader::skip_to(float target_t) {
  while (!eot()) {
    std::streampos frame_begin = fxtc.tellg();
    next_frame(true, true);

    if (t > target_t || std::abs(t - target_t) < teps) {
      fxtc.seekg(frame_begin);
      break;
    }
  }
}

void XTCReader::get_coordinates_3d(
    std::vector<std::array<float, 3>> &coords) const {
  coords.resize(nacur);
  for (uint32_t i = 0; i < nacur; i++) {
    coords[i][0] = X[3 * i];
    coords[i][1] = X[3 * i + 1];
    coords[i][2] = X[3 * i + 2];
  }
}
