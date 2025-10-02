/*
 * Example usage of the XTCReader C++ class
 */

#include "libxtc/xtc_reader.h"
#include "libxtc/xtc_writer.h"
#include <filesystem>
#include <iostream>

namespace fs = std::filesystem;

int main() {
  try {
    fs::path CURRENT = __FILE__;
    fs::path EXAMPLES_DIR = CURRENT.parent_path().parent_path() / "examples";
    fs::path md_file = EXAMPLES_DIR / "md.xtc";

    // Example 1: Basic usage - read all frames
    std::cout << "=== Example 1: Reading all frames ===" << std::endl;
    XTCReader reader(md_file);

    int frame_count = 0;
    while (!reader.eot()) {
      reader.next_frame();

      std::cout << "Frame " << frame_count << ", Time: " << reader.t << " ps"
                << ", Atoms: " << reader.get_natoms() << std::endl;

      // Access box vectors
      std::cout << "  Box: [" << reader.box[0][0] << ", " << reader.box[1][1]
                << ", " << reader.box[2][2] << "]" << std::endl;

      // Access coordinates (flat array)
      const float *coords = reader.get_coordinates();
      std::cout << "  First atom position: [" << coords[0] << ", " << coords[1]
                << ", " << coords[2] << "]" << std::endl;

      frame_count++;
    }

    std::cout << "\nTotal frames read: " << frame_count << "\n" << std::endl;

    // Example 2: Read with time window and sampling
    std::cout << "=== Example 2: Time window and sampling ===" << std::endl;
    XTCReader reader2(md_file,
                      100.0f, // start at 100 ps
                      500.0f, // end at 500 ps
                      10.0f); // read every 10 ps

    frame_count = 0;
    while (!reader2.eot()) {
      reader2.next_frame();
      std::cout << "Frame " << frame_count << ", Time: " << reader2.t << " ps"
                << std::endl;
      frame_count++;
    }

    std::cout << "\nFrames in time window: " << frame_count << "\n"
              << std::endl;

    // Example 3: Using 3D coordinate access
    std::cout << "=== Example 3: 3D coordinate access ===" << std::endl;
    XTCReader reader3(md_file);

    reader3.next_frame();

    std::vector<std::array<float, 3>> coords_3d;
    reader3.get_coordinates_3d(coords_3d);

    std::cout << "Number of atoms: " << coords_3d.size() << std::endl;
    std::cout << "First 3 atom positions:" << std::endl;
    for (size_t i = 0; i < std::min(size_t(3), coords_3d.size()); i++) {
      std::cout << "  Atom " << i << ": [" << coords_3d[i][0] << ", "
                << coords_3d[i][1] << ", " << coords_3d[i][2] << "]"
                << std::endl;
    }

    // Example 4: Calculate RMSD between frames
    std::cout << "\n=== Example 4: Computing with coordinates ===" << std::endl;
    XTCReader reader4(md_file);

    reader4.next_frame();
    std::vector<float> ref_coords = reader4.X; // Save first frame

    reader4.next_frame();
    const std::vector<float> &curr_coords = reader4.X;

    // Simple RMSD calculation (should align first, but this is just an example)
    double sum_sq = 0.0;
    for (size_t i = 0; i < ref_coords.size(); i++) {
      double diff = curr_coords[i] - ref_coords[i];
      sum_sq += diff * diff;
    }
    double rmsd = std::sqrt(sum_sq / (ref_coords.size() / 3.0));

    std::cout << "RMSD between frame 0 and frame 1: " << rmsd << " nm"
              << std::endl;

    // Example 5: Read only headers (fast scanning)
    std::cout << "\n=== Example 5: Fast header scanning ===" << std::endl;
    XTCReader reader5(md_file);

    frame_count = 0;
    float total_time = 0.0f;
    while (!reader5.eot()) {
      reader5.next_frame(true); // hdr_only = true
      total_time = reader5.t;
      frame_count++;
    }

    std::cout << "Scanned " << frame_count << " frames" << std::endl;
    std::cout << "Trajectory duration: " << total_time << " ps" << std::endl;

    // Example 6: Write XTC
    std::cout << "\n=== Example 6: Write XTC ===" << std::endl;
    // Create writer with precision 1000.0
    XTCWriter writer("output.xtc", 1000.0f);

    // Define simulation box (cubic, 5.0 nm)
    std::array<std::array<float, 3>, 3> box = {
        {{{5.0f, 0.0f, 0.0f}}, {{0.0f, 5.0f, 0.0f}}, {{0.0f, 0.0f, 5.0f}}}};

    // Number of atoms
    const uint32_t natoms = 1000;
    std::cout << "Writing " << natoms << " atoms" << std::endl;
    const int nframes = 100;
    std::cout << "Writing " << nframes << " frames" << std::endl;

    // Write 100 frames
    for (uint32_t frame = 0; frame < nframes; frame++) {
      float time = frame * 0.002f; // 2 ps timestep

      // Generate some example coordinates (random walk)
      std::vector<float> coords(natoms * 3);
      for (uint32_t i = 0; i < natoms; i++) {
        // Simple trajectory: atoms moving in a circle
        float angle = 2.0f * M_PI * i / natoms + time;
        float radius = 2.0f;
        coords[i * 3 + 0] = 2.5f + radius * std::cos(angle);
        coords[i * 3 + 1] = 2.5f + radius * std::sin(angle);
        coords[i * 3 + 2] = 2.5f + 0.1f * std::sin(time * 10.0f);
      }

      // Write frame with compression
      if (!writer.write_frame(frame, time, box, coords, true)) {
        std::cerr << "Failed to write frame " << frame << std::endl;
        return 1;
      }

      if (frame % 10 == 0) {
        std::cout << "Written frame " << frame << " at time " << time << " ps"
                  << std::endl;
      }
    }

    writer.close();
    std::cout << "Successfully written 100 frames to output.xtc" << std::endl;

    std::cout << "\n=== Example 7: Read written XTC ===" << std::endl;
    XTCReader new_reader("output.xtc");
    frame_count = 0;
    while (!new_reader.eot()) {
      new_reader.next_frame();
      frame_count++;
    }
    std::cout << "Read atoms: " << new_reader.get_natoms() << std::endl;
    std::cout << "Total frames: " << frame_count << std::endl;

  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
