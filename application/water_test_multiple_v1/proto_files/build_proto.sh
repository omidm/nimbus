#!/bin/bash

protoc --cpp_out=. ./physbam_vector_2d.proto
protoc --cpp_out=. ./physbam_range_2d.proto
protoc --cpp_out=. ./physbam_grid_2d.proto
protoc --cpp_out=. ./physbam_face_array_2d.proto
protoc --cpp_out=. ./app_face_array_2d.proto

protoc --cpp_out=. ./adv_vel_par.proto
