#!/usr/bin/python
import re
import sys
import os
import cgi
import time
sys.path.append("/n/curvature/data/pd")
from pd.common import CONNECT


######################################################################
# Write the HTML for a table
# headers cells have the optional form (colspan, width, header_text)
######################################################################
def HTML_Table(headers,rows):
    text=""
    text+=("<table width=\"100%\" cellspacing=0 cellpadding=2 border=0>\n")
    for h in headers:
        text+=("  <tr bgcolor=#7799ff>\n")
        for c in h:
            if type(c)==tuple:
                span,width,header=c
                text+=("    <th colspan=%d width=\"%s\">%s</th>\n"%(span,width,header))
            else:
                text+=("    <th>%s</th>\n"%c)
        text+=("  </tr>\n")
    color,alt_color="#eeeeee","#dddddd"
    for r in rows:
        color,alt_color=alt_color,color # alterante
        text+=("  <tr bgcolor=%s>\n"%color)
        for c in r:
            text+=("    <td>%s</td>\n"%c)
        text+=("  </tr>\n")
    
    text+=("</table>\n")
    return text



# grab fields
form=cgi.FieldStorage()
username_filter=""
state_filter=""

print "Content-type: text/html\n\n"

# get session list
client=CONNECT.sim_client()
sessions=client.Session_List()
client=None

# grab form options
if form.has_key("username_filter"): username_filter=form["username_filter"].value
if form.has_key("state_filter"): state_filter=form["state_filter"].value

# make the filter regular expressions
filter_conditions=[("username",re.compile(username_filter)),("state",re.compile(state_filter))]

# make a list of accepted sessions
accepted_sessions=[]
for i in sessions.keys():
    session=sessions[i]
    if reduce(lambda x,y: x and y,
              map(lambda x: x[1].match(session[x[0]]),filter_conditions)):
        accepted_sessions.append(i)
accepted_sessions.sort(lambda x,y: y-x)

# Build Python Table
user_fields_options={"last_frame": lambda x: str(x),
             "exit_time": lambda x: time.ctime(x),
             "exit_code": lambda x: str(x),
             "last_frame_duration": lambda x: "%.2f s"%x,
             "max_memory": lambda x: "%.2f MB"%x,
             "memory": lambda x: "%.2f MB"%x,
             "start_time": lambda x: time.ctime(x),
             "total_duration": lambda x: "%.2f s"%x
             }
used_user_fields=[]
for i in user_fields_options.keys():
   if form.has_key("field_"+i):
       used_user_fields.append(i)
user_fields=[(x,user_fields_options[x]) for x in used_user_fields]
labels=["ID","Username","State","Machine","Label"]
for field,format in user_fields:
    labels.append(field)
rows=[]
for i in accepted_sessions:
    session=sessions[i]
    row=[session["id"],session["username"],session["state"],session["machine"],session["label"]]
    for field,format in user_fields:
        if session["user_status"].has_key(field):
            row.append(format(session["user_status"][field]))
        else:
            row.append("N/A")
    #print repr(row)
    rows.append(row)

# Write page
print """<html><head>
<title>PD Simulation Status</title></head>
<link rel="stylesheet" type="text/css" href="style_status.css" />
<body>"""
print "<h1>PD Simulation Status</h1>"

print HTML_Table([labels],rows)
print "<p><form method=get action='WEB_STATUS.cgi'>"
print "State Filter: <input type=input name='state_filter' value='%s'><br>"%state_filter
print "Username Filter: <input type=input name='username_filter' value='%s'><br>"%username_filter
for i in user_fields_options.keys():
    checked=""
    if i in used_user_fields: checked="checked"
    print "<input type=checkbox name=\"field_%s\" %s> %s<br>"%(i,checked,i)
print "<input type=submit value='Update'><br>"
print "</form></p>"
print "</body></html>"
