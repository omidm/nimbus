#include "server.h"
#include <string>
#include <sstream>
#include <iostream>
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>

using boost::asio::ip::tcp;

Server::Server(unsigned int _connection_port_no, Scheduler* sch)
: scheduler(sch),
  connection_subscription_thread(NULL),
  connection_port_no(_connection_port_no)
{}

Server::~Server()
{
  boost::mutex::scoped_lock lock(map_mutex);
  for (ConnectionMapIter iter = connections.begin(); iter != connections.end();
      ++iter)
  {
    ServerConn* scon = iter->second;
    delete scon;
  }
  connections.clear();
}

void Server::receive_msg(const std::string& msg, ServerConn* conn)
{
  // FIXME: memory leak for killing and ending thread listening for new connections.

  std::cout<<"\nReceived msg "<<msg<<"\n";

  boost::mutex::scoped_lock lock(map_mutex);
  for (ConnectionMapIter iter = connections.begin(); iter != connections.end();
      ++iter)
  {
    if (iter->first != conn->get_id())
      iter->second->send_msg(msg);
  }
}

void Server::listen_for_new_connections()
{
  boost::asio::io_service io_service;
  tcp::acceptor acceptor(
      io_service, tcp::endpoint(tcp::v4(), connection_port_no));

  while (true)
  {
    tcp::socket* socket = new tcp::socket(io_service);
    acceptor.accept(*socket);
    {
      std::cout<<"\nCreating new connection\n";
      boost::mutex::scoped_lock lock(map_mutex);
      ServerConn* sc = new ServerConn(this,socket);
      connections[sc->get_id()] = sc;
      sc->start_listening();
    }
  }
}

void Server::run()
{
  connection_subscription_thread = new boost::thread(
      boost::bind(&Server::listen_for_new_connections,this));

  connection_subscription_thread->join();
}


ServerConn::ServerConn(Server* s, tcp::socket* sock) 
: server(s),
  socket(sock)
{
  static ConnectionId id_assigner = 0;
  id = id_assigner ++;
}


void ServerConn::listen_for_msgs()
{
  while (true)
  {
    boost::asio::streambuf response;
    boost::asio::read_until(*socket, response, ";");

    std::istream is(&response);
    std::string msg;
    std::getline(is, msg);
    server->receive_msg(msg, this);
  }
}

void ServerConn::start_listening()
{
  listening_thread = new boost::thread(
      boost::bind(&ServerConn::listen_for_msgs,this));
}

void ServerConn::send_msg(const std::string& msg)
{
  boost::system::error_code ignored_error;
  boost::asio::write(*socket, boost::asio::buffer(msg),
      boost::asio::transfer_all(), ignored_error);
}

ServerConn::~ServerConn()
{
  // FIXME: not actually cleaning up listening thread.
}



