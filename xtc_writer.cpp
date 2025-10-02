#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

class XTCWriter {
public:
  static constexpr int32_t MXTC = 1995;
  static constexpr float DEFAULT_PRECISION = 1000.0f; // 1/nm precision

  XTCWriter(const std::string &fname, int natoms,
            float precision = DEFAULT_PRECISION)
      : natoms_(natoms), precision_(precision), frame_index_(0) {
    file_.open(fname, std::ios::binary);
    if (!file_.is_open()) {
      throw std::runtime_error("Cannot create file: " + fname);
    }
  }

  ~XTCWriter() { close(); }

  // Write a frame
  void write_frame(float time, const std::vector<std::array<float, 3>> &coords,
                   const std::array<std::array<float, 3>, 3> &box) {
    if (static_cast<int>(coords.size()) != natoms_) {
      throw std::runtime_error("Coordinate array size does not match natoms");
    }

    write_be_int32(MXTC);
    write_be_uint32(natoms_);
    write_be_uint32(frame_index_);
    write_be_float(time);

    // Write box
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        write_be_float(box[i][j]);
      }
    }

    // For simplicity, write uncompressed for small systems
    // or always uncompressed (lsize <= 9 means uncompressed)
    bool write_uncompressed = (natoms_ <= 9);

    if (write_uncompressed) {
      write_be_uint32(natoms_); // lsize

      // Write coordinates
      for (int i = 0; i < natoms_; i++) {
        for (int j = 0; j < 3; j++) {
          write_be_float(coords[i][j]);
        }
      }
    } else {
      // Write compressed (requires pack_frame function from libxtc)
      // For now, we'll write uncompressed and note this limitation
      write_uncompressed_frame(coords);
    }

    frame_index_++;
  }

  // Simplified interface matching typical usage
  void write(float time, const std::vector<std::array<float, 3>> &coords,
             const std::array<std::array<float, 3>, 3> &box) {
    write_frame(time, coords, box);
  }

  // Alternative: write with 1D coordinate array
  void write_frame(float time,
                   const float *coords, // natoms * 3 array
                   const float *box) {  // 9 element array
    std::vector<std::array<float, 3>> coord_vec(natoms_);
    for (int i = 0; i < natoms_; i++) {
      coord_vec[i][0] = coords[i * 3 + 0];
      coord_vec[i][1] = coords[i * 3 + 1];
      coord_vec[i][2] = coords[i * 3 + 2];
    }

    std::array<std::array<float, 3>, 3> box_arr;
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        box_arr[i][j] = box[i * 3 + j];
      }
    }

    write_frame(time, coord_vec, box_arr);
  }

  int frame_count() const { return frame_index_; }

  void close() {
    if (file_.is_open()) {
      file_.close();
    }
  }

private:
  std::ofstream file_;
  int natoms_;
  float precision_;
  int frame_index_;

  void write_be_int32(int32_t val) {
    uint8_t buf[4];
    buf[0] = (val >> 24) & 0xFF;
    buf[1] = (val >> 16) & 0xFF;
    buf[2] = (val >> 8) & 0xFF;
    buf[3] = val & 0xFF;
    file_.write(reinterpret_cast<char *>(buf), 4);
  }

  void write_be_uint32(uint32_t val) {
    write_be_int32(static_cast<int32_t>(val));
  }

  void write_be_float(float val) {
    uint32_t i;
    std::memcpy(&i, &val, sizeof(float));
    write_be_uint32(i);
  }

  void
  write_uncompressed_frame(const std::vector<std::array<float, 3>> &coords) {
    write_be_uint32(natoms_); // lsize

    for (int i = 0; i < natoms_; i++) {
      for (int j = 0; j < 3; j++) {
        write_be_float(coords[i][j]);
      }
    }
  }
};
