#ifndef __SERVER_HPP__
#define __SERVER_HPP__

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <map>
#include "parser.h"

typedef unsigned int ConnectionId;
class Scheduler;
class ServerConn;

using boost::asio::ip::tcp;

class Server
{
public:
    Server(unsigned int _connection_port_no, Scheduler* sch);
    ~Server();

    Scheduler * scheduler;
    
    void run();
    void receive_msg(const std::string& msg, ServerConn* conn);

    
private:

    void listen_for_new_connections();
    
    typedef std::map<ConnectionId, ServerConn*> ConnectionMap;
    typedef ConnectionMap::iterator ConnectionMapIter;
    ConnectionMap connections;

    boost::mutex map_mutex;
    boost::thread* connection_subscription_thread;
    unsigned int connection_port_no;
};


class ServerConn
{
public:
    ServerConn(Server* s, tcp::socket* sock);
    ~ServerConn();
    void send_msg(const std::string& msg);
    void start_listening();
    ConnectionId get_id() const
    {
        return id;
    }
private:
    void listen_for_msgs();
    ConnectionId id;
    Server* server;
    tcp::socket* socket;
    boost::thread* listening_thread;
};

#endif
