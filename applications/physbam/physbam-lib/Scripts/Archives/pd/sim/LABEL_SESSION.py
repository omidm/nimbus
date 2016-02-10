#!/usr/bin/python
import sys
import os
import dialog
import CLIENT_LIBRARY

d=dialog.Dialog(dialog="dialog")
d.add_persistent_args(["--backtitle", "PSIM Label Session"])

id=CLIENT_LIBRARY.Get_Session(d)
CLIENT_LIBRARY.Label_Session(d,id)
