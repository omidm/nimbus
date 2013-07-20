#ifndef __SERVER_HPP__
#define __SERVER_HPP__

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <map>
#include "parser.h"

typedef unsigned int ConnectionId;
class Scheduler;
class SchedulerServerConnection;

using boost::asio::ip::tcp;

class SchedulerServer {
public:
    SchedulerServer(unsigned int _connection_port_no, Scheduler* sch);
    ~SchedulerServer();

    Scheduler * scheduler;
    
    void run();
    void receive_msg(const std::string& msg, SchedulerServerConnection* conn);

    
private:

    void listen_for_new_connections();
    
    typedef std::map<ConnectionId, SchedulerServerConnection*> ConnectionMap;
    typedef ConnectionMap::iterator ConnectionMapIter;
    ConnectionMap connections;

    boost::mutex map_mutex;
    boost::thread* connection_subscription_thread;
    unsigned int connection_port_no;
};


class SchedulerServerConnection {
public:
    SchedulerServerConnection(SchedulerServer* s, tcp::socket* sock);
    ~SchedulerServerConnection();
    void send_msg(const std::string& msg);
    void start_listening();
    ConnectionId get_id() const {
        return id;
    }
    
private:
    void listen_for_msgs();
    ConnectionId id;
    SchedulerServer* server;
    tcp::socket* socket;
    boost::thread* listening_thread;
};

#endif
