#!/usr/bin/python
import sys
import os
import dialog
import CLIENT_LIBRARY

d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle","PSIM Connect To Session"])

id=CLIENT_LIBRARY.Get_Session(d)

def Connect_Session(d,id):
    session_directory=CLIENT_LIBRARY.client.Session_Directory(id);
    info=CLIENT_LIBRARY.client.Session_Info(id)
    machine=info["machine"]
    if machine==None:
        print "Session does not have a machine associated with it. Activate it first or use GOTO_SANDBOX"
        sys.exit(0)
    sshcommand="ssh "+machine+" -t \"cd "+session_directory+"; export PSIM_SESSION_DIRECTORY="+session_directory+";export PSIM_SESSION_ID="+str(id)+";export PSIM_SERVER_PORT="+";bash\""
    print "Connecting To: "+machine+"  - Session ID: "+str(id)+"  - Label: "+info["label"]
    CLIENT_LIBRARY.client=None
    os.system(sshcommand);

Connect_Session(d,id);
