#!/usr/bin/python

# To effectively use this to change directories, use something like
# alias pdsgosim='`GET_SESSION_DIRECTORY.py 3>&1 1>&2 2>&3`'

import sys
import os
import dialog
import CLIENT_LIBRARY

d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle","PSIM Goto Sandbox"])

id=CLIENT_LIBRARY.Get_Session(d)

def Get_Session_Directory(d,id):
    session_directory=CLIENT_LIBRARY.client.Session_Directory(id)
    info=CLIENT_LIBRARY.client.Session_Info(id)
    sys.stderr.write("cd "+session_directory);

Get_Session_Directory(d,id)
sys.exit(0)
