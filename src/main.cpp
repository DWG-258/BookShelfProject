#include "file_reader.hpp"
#include <filesystem>
#include <iostream>
int main(int argc, char* argv[]) {
    std::cout << "Hello, World!" << std::endl;
  std::filesystem::path current_path = std::filesystem::current_path();
  std::string file("/home/ezio/ClassProject3/adaptec1");
  std::cout << "Reading file: " << file << std::endl;
  auto fr = std::make_shared<FileReader>(file);
  fr->print();
  std::cout << "File read successfully!" << std::endl;
  return 0;
}
