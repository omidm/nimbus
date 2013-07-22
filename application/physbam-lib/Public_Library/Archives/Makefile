default: PhysBAM PhysBAM_OpenGL

PhysBAM:
	$(MAKE) -f Makefile.PhysBAM

PhysBAM_OpenGL:
	$(MAKE) -f Makefile.PhysBAM_OpenGL

%:
	$(MAKE) -f Makefile.PhysBAM $@
	$(MAKE) -f Makefile.PhysBAM_OpenGL $@
