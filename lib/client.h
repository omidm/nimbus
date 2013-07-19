#ifndef __SERVER_HPP__
#define __SERVER_HPP__

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <map>

typedef uint ConnectionId;

using boost::asio::ip::tcp;

class Client
{
public:
    Client(uint _connection_port_no);
    ~Client();
    
    void run();
    
    void receive_msg();

    void send_msg();

    
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



#endif
