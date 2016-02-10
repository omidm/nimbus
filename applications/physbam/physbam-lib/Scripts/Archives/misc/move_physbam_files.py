#!/usr/bin/python
import os
import sys
import utilities.file_utilities
import re
from optparse import OptionParser

def run_command(command_string,verbosity=1,noop=False):
    if(verbosity):
        print "running: ",command_string

    if(not noop): os.system(command_string)

def move_file(old_file_relative_path,new_file_relative_path,library_path,rename,noop=False):
    print "moving file %s to %s"%(old_file_relative_path,new_file_relative_path)
    old_class_name,extension=os.path.basename(old_file_relative_path).rsplit('.',1)
    new_class_name,extension=os.path.basename(new_file_relative_path).rsplit('.',1)
    old_file_abs_path=os.path.abspath(old_file_relative_path)
    new_file_abs_path=os.path.abspath(new_file_relative_path)
    library_path=os.path.abspath(library_path)
    library_path=library_path+"/"
    projects_path=library_path+"../Projects/"
    tools_path=library_path+"../Tools/"

    old_file_library_relative_path=old_file_abs_path.replace(library_path,'',1)
    new_file_library_relative_path=new_file_abs_path.replace(library_path,'',1)
    new_directory_library_relative_path=os.path.dirname(new_file_library_relative_path)

    print "\n\n*************Arguments*********************\n\n"
    print "old_file_abs_path=",old_file_abs_path
    print "new_file_abs_path=",new_file_abs_path

    print "old_file_library_relative_path=",old_file_library_relative_path
    print "new_file_library_relative_path=",new_file_library_relative_path

    print "old_class_name=",old_class_name
    print "new_class_name=",new_class_name
    print "\n\n*******************************************\n\n"

    #command_cvs_add_directory="""cvs add %s"""%new_directory_library_relative_path
    command_copy_file="""cp %s %s"""%(old_file_library_relative_path,new_file_library_relative_path)
    command_cvs_remove_old_file="""cvs remove -f %s"""%old_file_library_relative_path
    command_cvs_add_new_file="""cvs add %s"""%new_file_library_relative_path
    grep_dot_h_and_cpp="find ./ -iname *.h -o -iname *.cpp|xargs grep "
    command_search_replace_includes_to_new_path="""sed -i 's,%s,%s,g' `%s -rl %s`"""%(old_file_library_relative_path,new_file_library_relative_path,grep_dot_h_and_cpp,old_file_library_relative_path)
    command_search_replace_class_name="""sed -i 's,\<%s\>,%s,g' `%s -rlw %s`"""%(old_class_name,new_class_name,grep_dot_h_and_cpp,old_class_name)
    command_search_replace_file_ifndef="""sed -i 's,\<__%s__\>,__%s__,g' `%s -rl __%s__`"""%(old_class_name,new_class_name,grep_dot_h_and_cpp,old_class_name)

    current_working_directory=os.getcwd()

    print "\n\nchanging directory to %s"%library_path
    os.chdir(library_path)
    #if(not noop): utilities.file_utilities.make_directory(new_directory_library_relative_path)
    #run_command(command_cvs_add_directory,noop=noop)
    run_command(command_copy_file,noop=noop)
    run_command(command_cvs_remove_old_file,noop=noop)
    run_command(command_cvs_add_new_file,noop=noop)
    run_command(command_search_replace_includes_to_new_path,noop=noop)
    if(rename): run_command(command_search_replace_class_name,noop=noop)
    if(rename): run_command(command_search_replace_file_ifndef,noop=noop)

    print "\n\nchanging directory to %s"%projects_path
    os.chdir(projects_path)
    run_command(command_search_replace_includes_to_new_path,noop=noop)
    if(rename): run_command(command_search_replace_class_name,noop=noop)

    print "\n\nchanging directory to %s"%tools_path
    os.chdir(tools_path)
    run_command(command_search_replace_includes_to_new_path,noop=noop)
    if(rename): run_command(command_search_replace_class_name,noop=noop)

    print "\n\nchanging directory to %s"%current_working_directory
    os.chdir(current_working_directory)

def main():
    usage="\n%prog library_path new_directory_relative_path <old_files_relative_path_list>\nor\n%prog --rename library_path new_file_relative_path old_file_relative_path"
    parser=OptionParser(usage)
    parser.add_option("-r","--rename",action="store_true",dest="rename",help="file is being renamed")
    parser.add_option("-n","--noop",action="store_true",dest="noop",help="only print commands to be run. Don't actually run the commands")

    (options,args)=parser.parse_args(sys.argv[1:])

    if(options.rename):
        print "renaming,len=",len(args)
        if(len(args)!=3): parser.error("incorrect number of arguments")
    else:
        if(len(args)<3): parser.error("incorrect number of arguments")

    library_path=args[0]
    print "library_path=",library_path
    if(options.rename):
        new_file_relative_path=args[1]
        old_file_relative_path=args[2]
        move_file(old_file_relative_path,new_file_relative_path,library_path,rename=True,noop=options.noop)
    else:
        new_directory_relative_path=args[1]
        filename_regex=re.compile(".*\.(h|cpp)$")
        match=filename_regex.match(new_directory_relative_path)
        if(match):
            print "new_directory_relative_path=%s can't be a directory"%new_directory_relative_path
            sys.exit(0)

        for old_file_relative_path in args[2:]:
            new_file_relative_path=os.path.join(new_directory_relative_path,os.path.basename(old_file_relative_path))
            move_file(old_file_relative_path,new_file_relative_path,library_path,rename=False,noop=options.noop)
            print "\n\n"

if __name__=="__main__":
    main()



# e.g.
#move_physbam_files.py /data/kwatra/PhysBAM/Public_Library/ PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Rigids_Components/ PhysBAM_Solids/PhysBAM_Rigids/OpenGL_Components/*.{h,cpp}
#move_physbam_files.py /data/kwatra/PhysBAM/Public_Library/ PhysBAM_Rendering/PhysBAM_OpenGL_Solids/OpenGL_Deformable_Components/ PhysBAM_Solids/PhysBAM_Deformables/OpenGL_Components/*.{h,cpp}
#move_physbam_files.py /data/kwatra/PhysBAM/Public_Library/ PhysBAM_Rendering/PhysBAM_OpenGL_Fluids/OpenGL_Incompressible_Components/ PhysBAM_Fluids/PhysBAM_Incompressible/OpenGL_Components/*.{h,cpp}
#move_physbam_files.py /data/kwatra/PhysBAM/Public_Library/ PhysBAM_Rendering/PhysBAM_OpenGL_Dynamics/OpenGL_Dynamics_Components/ PhysBAM_Dynamics/OpenGL_Components/*.{h,cpp}
