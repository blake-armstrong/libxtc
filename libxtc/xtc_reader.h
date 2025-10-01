/*
 * Copyright (c) 2020, Nikolay A. Krylov
 * All rights reserved.
 *
 * C++ XTCReader class - header file
 */

#ifndef XTC_READER_H
#define XTC_READER_H

#include "xtc_unpack.h"
#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

class XTCReader {
public:
  // Main data arrays
  std::vector<float> X;                    // Atom coordinates [natoms * 3]
  std::array<std::array<float, 3>, 3> box; // Unit cell vectors [3x3]
  float t;                                 // Current frame time

  // Constructor
  XTCReader(const std::string &fname, float beg = -1.0f, float end = -1.0f,
            float dt = -1.0f, int nt = 4);

  // Destructor
  ~XTCReader();

  // Read next frame from trajectory
  // hdr_only: if true, only read header and skip coordinates
  // no_check: if true, skip time checking for dt sampling
  bool next_frame(bool hdr_only = false, bool no_check = false);

  // Check if end of trajectory reached
  bool eot();

  // Get current number of atoms
  uint32_t get_natoms() const { return nacur; }

  // Get coordinates as pointer (for C-style access)
  const float *get_coordinates() const { return X.data(); }

  // Get coordinates reshaped as Nx3 (utility method)
  void get_coordinates_3d(std::vector<std::array<float, 3>> &coords) const;

private:
  // Constants
  static constexpr int32_t mxtc = 1995;
  static constexpr int hdr1 = 4 * 4 + 4 * 9 + 4;
  static constexpr int hdr2 = 4 * 9;
  static constexpr int hdrfull = hdr1 + hdr2;
  static constexpr int flt3len = 3 * 4;
  static constexpr int alignmentBytes = 4;
  static constexpr int alignmentBytesMinusOne = alignmentBytes - 1;

  // File handling
  std::ifstream fxtc;
  std::streampos fsize;
  std::string filename;

  // Frame data
  frame_data frd;
  std::vector<uint8_t> buf;
  uint32_t nacur;
  uint32_t lsize;
  int datalen;

  // Trajectory parameters
  float beg, end;
  float dt;
  int iframe;
  float t0;
  float dtrrj; // Time between consecutive frames in file
  float teps;  // Time epsilon for comparisons
  bool bigdt;  // Flag for dt sampling strategy

  // Helper methods
  static int align4(int l) {
    return (l + alignmentBytesMinusOne) & ~alignmentBytesMinusOne;
  }

  bool read_header(bool check = false);
  void decompress();
  bool search_frame_hdr();
  void search_t(float target_t);
  void skip_to(float target_t);
  float tnext() const;

  // XDR unpacking utilities
  int32_t unpack_int(const uint8_t *&ptr);
  uint32_t unpack_uint(const uint8_t *&ptr);
  float unpack_float(const uint8_t *&ptr);

  // Endian conversion
  template <typename T> T byteswap(T value);
};

#endif // XTC_READER_H
