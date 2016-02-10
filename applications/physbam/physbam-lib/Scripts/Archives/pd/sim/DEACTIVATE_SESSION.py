#!/usr/bin/python
import sys
import os
import dialog
import CLIENT_LIBRARY

d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle", "PSIM Deactivate Session"])

id=CLIENT_LIBRARY.Get_Session(d,query_states=["active"])
CLIENT_LIBRARY.Deactivate_Session(d,id)
