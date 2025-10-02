/*
 * Copyright (c) 2025
 * All rights reserved.
 *
 * C++ XTCWriter class - header file
 */

#ifndef XTC_WRITER_H
#define XTC_WRITER_H

#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

class XTCWriter {
public:
  // Constructor
  // fname: output filename
  // precision: compression precision (default 1000.0)
  XTCWriter(const std::string &fname, float precision = 1000.0f);

  // Destructor
  ~XTCWriter();

  // Write a frame to the trajectory
  // natoms: number of atoms
  // step: frame/step number
  // time: simulation time
  // box: 3x3 box matrix
  // coords: atom coordinates [natoms * 3]
  // compress: whether to use compression (default true)
  bool write_frame(uint32_t natoms, uint32_t step, float time,
                   const std::array<std::array<float, 3>, 3> &box,
                   const float *coords, bool compress = true);

  // Convenience method with vector input
  bool write_frame(uint32_t step, float time,
                   const std::array<std::array<float, 3>, 3> &box,
                   const std::vector<float> &coords, bool compress = true);

  // Convenience method with 3D coordinate input
  bool write_frame(uint32_t step, float time,
                   const std::array<std::array<float, 3>, 3> &box,
                   const std::vector<std::array<float, 3>> &coords,
                   bool compress = true);

  // Flush and close the file
  void close();

  // Check if file is open
  bool is_open() const { return fxtc.is_open(); }

private:
  // Constants
  static constexpr int32_t mxtc = 1995;
  static constexpr int alignmentBytes = 4;
  static constexpr int alignmentBytesMinusOne = alignmentBytes - 1;

  // File handling
  std::ofstream fxtc;
  std::string filename;

  // Compression parameters
  float precision;

  // Temporary buffers
  std::vector<uint8_t> compressed_buffer;
  std::vector<int32_t> int_coords;

  // Helper methods
  static int align4(int l) {
    return (l + alignmentBytesMinusOne) & ~alignmentBytesMinusOne;
  }

  // XDR packing utilities
  void pack_int(std::vector<uint8_t> &buf, int32_t value);
  void pack_uint(std::vector<uint8_t> &buf, uint32_t value);
  void pack_float(std::vector<uint8_t> &buf, float value);

  // Compression methods
  bool compress_coords(uint32_t natoms, const float *coords,
                       std::vector<uint8_t> &output, int32_t &nbytes);

  // Endian conversion
  template <typename T> T byteswap(T value);

  // Calculate integer coordinates and bounds
  void calculate_bounds(uint32_t natoms, const float *coords, int32_t minint[3],
                        int32_t maxint[3], std::vector<int32_t> &int_coords);
};

#endif // XTC_WRITER_H
