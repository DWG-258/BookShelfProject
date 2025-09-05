#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include "placedata.hpp"

class FileReader{
    public:
        /**
         * @brief Construct a new File Reader object
         * 
         * @param filepath 需要读取的文件夹
         * @return none
         * 
         */
        FileReader(const std::string& filepath){
            
        }

    private:
    void read_aux();
    void read_nets();
    void read_nodes();
    void read_pl();
    void read_scl();
    void read_wts();
};