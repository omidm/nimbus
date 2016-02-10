from pd.common import CONFIG
from pd.common import SOCKET
import sys
import time
import os
import socket
import dialog

client=None
try:
    client=SOCKET.CLIENT(CONFIG.pdsim_server_host,CONFIG.pdsim_server_port,
                  (CONFIG.client_private_key_file,CONFIG.client_certificate_file,CONFIG.ca_certificate_file))

#    client=SOCKET.CLIENT(os.environ["PSIM_SERVER_HOST"],int(os.environ["PSIM_SERVER_PORT"]),
#                         (os.environ["PSIM_CLIENT_KEY"],os.environ["PSIM_CLIENT_CERT"],os.environ["PSIM_CA_CERT"]))
except KeyError:
    print "You must define the environment variables PSIM_SERVER_HOST and PSIM_SERVER_PORT, PSIM_CLIENT_KEY, \nPSIM_CA_CERT, PSIM_CLIENT_CERT"
    sys.exit(1)
except socket.error:
    print "Unable to connect to server"
    sys.exit(1)

username=os.environ["USER"]


def handle_dialog_code(d,code):
    if code in (d.DIALOG_CANCEL, d.DIALOG_ESC):
        sys.exit(1)
    return 0


def New_Session(d):

    # Get memory
    (code,answer)=d.inputbox("How much memory (GBytes) will this job take?",init="4")
    handle_dialog_code(d,code)
    memory=int(answer)

    # Get CPUs
    (code,answer)=d.inputbox("How many CPUs will this job take?",init="1")
    handle_dialog_code(d,code)
    cpus=int(answer)

    # Get Label
    (code,label)=d.inputbox("What do you call this job?",init="")
    handle_dialog_code(d,code)

    # Try making the job
    session_info=client.Create_Session(username,memory,cpus)
    id=None
    try:
        id=session_info["id"]
    except:
        d.infobox("Failed to get session")
        sys.exit(1)
    client.Label_Session(id,label)

    return id

# query_states is a list of states we want
# users is a list of users we want
def Get_Session(d,query_states=None,users=None):
    sessions=client.Session_List()
    formatted_choices=[]
    for i in sessions.keys():
        session=sessions[i]
        if query_states and not session["state"] in query_states: continue
        if users and not session["username"] in users: continue
        statestr=session["state"]
        if statestr=="active": statestr="active@%s"%session["machine"]
        formatted_choices.append((str(session["id"]),"%-8s %-20s %-200s"%(session["username"],statestr,session["label"])))
    if len(formatted_choices)==0:
        print "No sessions found with state %s and users %s"%(repr(query_states),repr(users))
        sys.exit(1)
    (code,session)=d.menu("Choose session", width=230,height=-1,menu_height=0,choices=formatted_choices)
    handle_dialog_code(d,code)
    return int(session)

def Label_Session(d,id):
    label=client.Session_Info(id)["label"]
    (code,label)=d.inputbox("What should the label be?",init=label)
    handle_dialog_code(d,code)
    client.Label_Session(id,label)

def Deactivate_Session(d,id):
    (code,state)=d.menu("What state should it get?",width=60,height=-1,menu_height=0,choices=[("inactive","Not running but might be soon"),("done","Pretty much done")])
    handle_dialog_code(d,code)
    client.Deactivate_Session(id,state)

def Activate_Session(d,id):
    while 1:
        hosts=client.Host_List()
        formatted=[]
        for host in hosts.keys():
            claims=hosts[host]["claims"]
            mem=hosts[host]["max_memory"]
            cpu=hosts[host]["max_cpus"]
            users=[]
            for claim in claims.keys():
                cpu-=claims[claim]["cpus"]
                mem-=claims[claim]["memory"]
                users.append(claims[claim]["user"])
            avail_string="Free CPU=%2d Free Mem=%3d claims=%s"%(cpu,mem,", ".join(users))
            formatted.append((host,avail_string))
        if len(formatted)==0:
            print "No hosts"
            sys.exit(1)
        (code,hostname)=d.menu("Which host you would like?", width=60,height=-1,menu_height=0,choices=formatted)
        handle_dialog_code(d,code)

        # now try to attach
        try:
            client.Activate_Session(id,hostname)
            d.infobox("Session %d successfully attached to %s"%(id,hostname));
            break
        except SOCKET.COMMAND_EXCEPTION,e:
            d.msgbox("The following error occured when trying to get host:\n\n"+str(e))

def Status(d,id):
    info=client.Session_Info(id)
    lines=["%20s %d"%("ID:",info["id"]),
           "%20s %s"%("Label:",info["label"]),
           "%20s %s"%("State:",info["state"]),
           "%20s %d"%("Memory (GB):",info["memory"]),
           "%20s %d"%("CPUs:",info["cpus"]),
           "%20s %s"%("Username:",info["username"]),
           "%20s %s"%("Machine:",info["machine"]),
           "%20s %s"%("Date Created:",time.ctime(info["created_date"])),
           "",
           "User Status",
           "-----------"]
    status=info["user_status"]
    stats=status.keys()
    stats.sort()
    for i in stats:
        lines.append("   %20s : %s"%(i,repr(status[i])))

    d.msgbox(width=-1,height=-1,text="\n".join(lines))
    
    
    
