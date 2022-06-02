#include <iostream>
#include <algorithm>
#include "SeqIO.hh"

int main(int argc, char** argv){
    SeqIO::Reader<> in(argv[1]);
    SeqIO::Writer<> out(argv[2]);

    while(true){
        int64_t len = in.get_next_read_to_buffer();
        if(len == 0) break;
        std::reverse(in.read_buf, in.read_buf + len);
        for(int64_t i = 0; i < len; i++) in.read_buf[i] = SeqIO::get_rc(in.read_buf[i]);
        out.write_sequence(in.read_buf, len);
    }
}