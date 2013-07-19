#include "client.h"
#include <string>
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <sstream>
#include <iostream>

using boost::asio::ip::tcp;

Client::Client(uint _connection_port_no)
: connection_port_no(_connection_port_no)
{
  io_service = new boost::asio::io_service();
  socket = new tcp::socket(*io_service);
}

Client::~Client()
{}

void Client::receive_msg()
{
  while (true)
  {
    boost::asio::streambuf response;
    boost::asio::read_until(*socket, response, ";");

    std::istream is(&response);
    std::string msg;
    std::getline(is, msg);
    std::cout<<"\nReceived msg: "<<msg<<"\n";
  }
}

void Client::send_msg()
{
  while(true)
  {
    std::string msg;
    std::getline(std::cin, msg);
    boost::system::error_code ignored_error;
    boost::asio::write(*socket, boost::asio::buffer(msg),
        boost::asio::transfer_all(), ignored_error);
  }
}

// Note: run from connection_subscription_thread
void Client::create_new_connections()
{
  tcp::resolver resolver(*io_service);
  tcp::resolver::query query("127.0.0.1", boost::to_string(connection_port_no));
  tcp::resolver::iterator iterator = resolver.resolve(query);
  boost::system::error_code error = boost::asio::error::host_not_found;
  socket->connect(*iterator,error);
}

void Client::run()
{
  create_new_connections();
  sending_thread = new boost::thread(boost::bind(&Client::send_msg,this));
  receiving_thread = new boost::thread(boost::bind(&Client::receive_msg,this));

  // for now, have no other work to do, so just wait until the listening
  // thread terminates.
  sending_thread->join();
  receiving_thread->join();
}


