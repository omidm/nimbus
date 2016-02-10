This is where libraries that are used in PhysBAM should go. Ideally, 
one should place the sources for them here and integrate them 
into the CMake build system as an ExternalProject so that they 
are built automatically when needed.

If you need to edit library-specific makefiles, use 

"PhysBAM Edit: <description>"

as a comment in the appropriate file.