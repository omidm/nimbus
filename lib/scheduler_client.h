#ifndef __SCHEDULER_CLIENT_HPP__
#define __SCHEDULER_CLIENT_HPP__

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <map>
#include "lib/scheduler_command.h"

typedef uint ConnectionId;

using boost::asio::ip::tcp;

class SchedulerClient {
public:
    SchedulerClient(uint _connection_port_no);
    ~SchedulerClient();
    
    void run();
    SchedulerCommand* receiveCommand();
    void sendCommand(SchedulerCommand* c);

private:
    void create_new_connections();
    
    boost::asio::io_service* io_service;
    
    //socket for connection
    tcp::socket* socket;

    //thread for receiving messages
    boost::thread* receiving_thread;

    //thread for sending messages
    boost::thread* sending_thread;

    // port number to connect to the server
    ConnectionId connection_port_no;
};



#endif  // __SCHEDULER_CLIENT_HPP__
