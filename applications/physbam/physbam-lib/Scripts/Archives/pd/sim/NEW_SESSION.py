#!/usr/bin/python
import sys
import os
import dialog
import CLIENT_LIBRARY

d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle", "pythondialog demo"])

id=CLIENT_LIBRARY.New_Session(d)
code=d.yesno("Job created successfully!\n Activate Session?")
if not code:
    CLIENT_LIBRARY.Activate_Session(d,id)
