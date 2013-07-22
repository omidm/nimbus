FUNCTION(Create_Source_Group groupName fullPathToGroupParent sourceFiles)
	FOREACH(currentFile ${ARGV2};${ARGN})
		FILE(RELATIVE_PATH folder ${fullPathToGroupParent}/${groupName} ${currentFile})
		GET_FILENAME_COMPONENT(filename ${folder} NAME)
		IF(NOT folder STREQUAL "")
			STRING(REPLACE "${filename}" "" folderlast ${folder})
			IF(NOT folderlast STREQUAL "")
				STRING(REGEX REPLACE "/+$" "" folderlast ${folderlast})
				IF(NOT folderlast STREQUAL "")
					STRING(REPLACE "/" "\\" folderlast ${folderlast})
				ENDIF(NOT folderlast STREQUAL "")
			ENDIF(NOT folderlast STREQUAL "")
			SOURCE_GROUP("\\${folderlast}" FILES ${currentFile})
		ENDIF(NOT folder STREQUAL "")	
	ENDFOREACH(currentFile ${ARGV2};${ARGN})
ENDFUNCTION(Create_Source_Group)