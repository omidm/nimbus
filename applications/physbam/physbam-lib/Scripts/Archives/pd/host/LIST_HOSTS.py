#!/usr/bin/python


from pd.common import CONNECT

client=CONNECT.host_client()

hosts=client.Host_List().items()
print "/-----------------+------------\\"
print "| %-15s | %-10s |"%("Host","User")
print "|-----------------+------------|"
hosts.sort(lambda x,y:cmp(x[0],y[0]))
for host in hosts:
   hostname=host[0]
   user=""
   if host[1].has_key("user"):
      user=host[1]["user"]
      if user==None:
        user=""
   print "| %-15s | %-10s |"%(hostname,user)
print "\\-----------------+------------/"


