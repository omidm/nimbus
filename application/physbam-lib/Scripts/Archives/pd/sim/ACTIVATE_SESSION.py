#!/usr/bin/python
import sys
import os
import dialog
import CLIENT_LIBRARY


d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle", "PSIM Attach Session"])

id=CLIENT_LIBRARY.Get_Session(d,query_states=["inactive","done"])
CLIENT_LIBRARY.Activate_Session(d,id)
