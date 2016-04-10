/*
 * Copyright 2013 Stanford University.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * - Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * - Neither the name of the copyright holders nor the names of
 *   its contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
 * THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

 /*
  * Server side of the Nimbus scheduler protocol.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#include "src/shared/scheduler_server.h"
#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>
#include <boost/tokenizer.hpp>
#include <string>
#include <sstream>
#include <iostream>  // NOLINT
#include "src/shared/dbg.h"

#define SERVER_TCP_SEND_BUF_SIZE 10485760  // 10MB
#define SERVER_TCP_RECEIVE_BUF_SIZE 10485760  // 10MB

using boost::asio::ip::tcp;

namespace nimbus {

SchedulerServer::SchedulerServer(port_t port_no)
  : listening_port_(port_no) {
    bouncer_thread_active_ = false;
    total_bytes_sent_ = 0;
    total_bytes_received_ = 0;
  }

SchedulerServer::~SchedulerServer() {
  {
    boost::mutex::scoped_lock lock(worker_mutex_);
    for (SchedulerWorkerList::iterator iter = workers_.begin();
         iter != workers_.end();
         ++iter)   {
      SchedulerWorker* worker = *iter;
      worker->MarkDead();
      delete worker;
    }
    workers_.clear();
  }
  delete acceptor_;
  delete io_service_;
}

bool SchedulerServer::Initialize() {
  io_service_ = new boost::asio::io_service();
  acceptor_ = new tcp::acceptor(*io_service_,
      tcp::endpoint(tcp::v4(), listening_port_));
  return true;
}

bool SchedulerServer::ReceiveCommands(SchedulerCommandList* storage,
                                      size_t maxCommands) {
  storage->clear();
  boost::mutex::scoped_lock lock(command_queue_mutex_);
  uint32_t pending = received_commands_.size();
  if (pending == 0) {
    return false;
  } else if (pending < maxCommands) {
    maxCommands = pending;
  }
  for (uint32_t i = 0; i < maxCommands; i++) {
    SchedulerCommand* command = received_commands_.front();
    received_commands_.pop_front();
    dbg(DBG_NET, "Copying command %s to user buffer.\n", command->ToString().c_str());
    storage->push_back(command);
  }
  return true;
}


bool SchedulerServer::ReceiveJobDoneCommands(JobDoneCommandList* storage,
                                             size_t maxCommands) {
  storage->clear();
  boost::mutex::scoped_lock lock(command_queue_mutex_);
  uint32_t pending = received_job_done_commands_.size();
  if (pending == 0) {
    return false;
  } else if (pending < maxCommands) {
    maxCommands = pending;
  }
  for (uint32_t i = 0; i < maxCommands; i++) {
    JobDoneCommand* command = received_job_done_commands_.front();
    received_job_done_commands_.pop_front();
    dbg(DBG_NET, "Copying job done command %s to user buffer.\n", command->ToString().c_str());
    storage->push_back(command);
  }
  return true;
}


void SchedulerServer::SendCommand(SchedulerWorker* worker,
                                  SchedulerCommand* command) {
  dbg(DBG_NET, "Sending command %s .\n", command->ToString().c_str());
  std::string data = command->ToNetworkData();
  SchedulerCommand::length_field_t len;
  len = htonl((uint32_t)data.length() + sizeof(len));
  std::string msg;
  msg.append((const char*)&len, sizeof(len));
  msg.append(data.c_str(), data.length());

  total_bytes_sent_ += msg.size();

  boost::mutex::scoped_lock lock(*(worker->connection()->mutex()));
  SchedulerServerConnection* connection = worker->connection();
  boost::system::error_code ignored_error;
  // Why are we IGNORING ERRORS!??!?!?
  boost::asio::write(*(connection->socket()), boost::asio::buffer(msg),
                     boost::asio::transfer_all(), ignored_error);
}

void SchedulerServer::SendCommands(SchedulerWorker* worker,
                                    SchedulerCommandList* commands) {
  SchedulerCommandList::iterator iter = commands->begin();
  for (; iter != commands->end(); ++iter) {
    SchedulerCommand* command = *iter;
    SendCommand(worker, command);
  }
}

void SchedulerServer::BroadcastCommand(SchedulerCommand* command) {
  SchedulerWorkerList::iterator it = workers_.begin();
  for (; it != workers_.end(); ++it) {
    SchedulerWorker* w = *it;
    SendCommand(w, command);
  }
}

void SchedulerServer::BroadcastCommands(SchedulerCommandList* commands) {
  SchedulerWorkerList::iterator it = workers_.begin();
  for (; it != workers_.end(); ++it) {
    SchedulerWorker* w = *it;
    SendCommands(w, commands);
  }
}

void SchedulerServer::ListenForNewConnections() {
  dbg(DBG_NET, "Scheduler server listening for new connections.\n");
  tcp::socket* socket = new tcp::socket(*io_service_);
  SchedulerServerConnection* server_connection =
    new SchedulerServerConnection(socket);
  acceptor_->async_accept(*server_connection->socket(),
                          boost::bind(&SchedulerServer::HandleAccept,
                                      this,
                                      server_connection,
                                      boost::asio::placeholders::error));
}

void SchedulerServer::HandleAccept(SchedulerServerConnection* connection,
                                   const boost::system::error_code& error) {
  if (!error) {
    dbg(DBG_NET, "Scheduler accepted new connection.\n");
    SchedulerWorker* worker =  AddWorker(connection);
    // Set the tcp send and receive buf size.
    // Note: you may have to increase the OS limits first.
    // Look at the nimbus/scripts/configure_tcp.sh for help.
    boost::asio::socket_base::send_buffer_size s_option(SERVER_TCP_SEND_BUF_SIZE);
    boost::asio::socket_base::receive_buffer_size r_option(SERVER_TCP_RECEIVE_BUF_SIZE);
    connection->socket()->set_option(s_option);
    connection->socket()->set_option(r_option);
    // Turn of Nagle algorithm.
    boost::asio::ip::tcp::no_delay nd_option(TCP_NODELAY_OPTION);
    connection->socket()->set_option(nd_option);
    ListenForNewConnections();
    boost::asio::async_read(*(worker->connection()->socket()),
                            boost::asio::buffer(worker->connection()->read_buffer(),
                                                worker->connection()->max_read_buffer_length()),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerServer::HandleRead,
                                        this,
                                        worker,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  } else {
    delete connection;
  }
}

using boost::tokenizer;
using boost::char_separator;

  /** Reads commands from buffer with size bytes. Puts commands into
   * internal command list for later reads. Returns the number of
   * bytes read into commands. If this is less than size, the
   * caller must save the unread bytes. */
size_t SchedulerServer::EnqueueCommands(char* buffer, size_t size) {
  size_t total_read = 0;
  char* read_ptr = buffer;
  size_t bytes_remaining = size;

  dbg(DBG_NET, "Read %i bytes from worker.\n", size);
  size_t header_len = sizeof(SchedulerCommand::length_field_t);

  while (bytes_remaining >= header_len) {
    // Read the header length field
    SchedulerCommand::length_field_t len;
    char* ptr = reinterpret_cast<char*>(&len);
    memcpy(ptr, read_ptr, header_len);
    len = ntohl(len);

    if (bytes_remaining < len) {  // Don't have a complete command
      return total_read;
    } else {
      std::string input(read_ptr + header_len, len - header_len);
      SchedulerCommand* command;
      if (SchedulerCommand::GenerateSchedulerCommandChild(input,
                                                          worker_command_table_,
                                                          command)) {
        dbg(DBG_NET, "Enqueueing command %s.\n", command->ToString().c_str());

        boost::mutex::scoped_lock lock(command_queue_mutex_);
        received_commands_.push_back(command);
        if ((command->type() == SchedulerCommand::JOB_DONE) &&
            bouncer_thread_active_) {
          JobDoneCommand *comm = reinterpret_cast<JobDoneCommand*>(command);
          JobDoneCommand* dup_comm = new JobDoneCommand(comm->job_id(),
                                                        comm->run_time(),
                                                        comm->wait_time(),
                                                        comm->max_alloc(),
                                                        comm->final(),
                                                        comm->mark_stat());
          received_job_done_commands_.push_back(dup_comm);
        }
      } else {
        dbg(DBG_NET, "Ignored unknown command: %s.\n", input.c_str());
        exit(-1);
      }
      dbg(DBG_NET, "Read %i bytes of %i available.\n", len, bytes_remaining);
      bytes_remaining -= len;
      total_read += len;
      read_ptr += len;
    }
  }

  return total_read;
}

  /* 
  for (size_t i = 0; i < size; i++) {
    if (buffer[i] == ';') {
      // There could be null character (\0) in the string before semicolon.
      buffer[i] = '\0';
      std::string input(start_pointer, string_len);

      SchedulerCommand* command;
      if (SchedulerCommand::GenerateSchedulerCommandChild(
          input, worker_command_table_, command)) {
        dbg(DBG_NET, "Adding command %s to queue.\n", command->ToNetworkData().c_str());
        boost::mutex::scoped_lock lock(command_queue_mutex_);
        received_commands_.push_back(command);
      } else {
        dbg(DBG_NET, "Ignored unknown command: %s.\n", input.c_str());
        exit(-1);
      }
      // Next string starts after the semicolon
      start_pointer = buffer + i + 1;
      string_len = 0;
    } else {
      ++string_len;
    }
  }
  // We've read this many bytes successfully into commands
  return start_pointer - buffer;
}

  */
void SchedulerServer::HandleRead(SchedulerWorker* worker,
                                 const boost::system::error_code& error,
                                 size_t bytes_transferred) {
  if (error) {
    dbg(DBG_NET|DBG_ERROR,
        "\n*\n*\n****************** Error %s receiving %i bytes from worker %i.\n*\n*\n",
        error.message().c_str(), bytes_transferred, worker->worker_id());

    // Signal controller about the down worker bu pushing a nitification.
    SchedulerCommand *command = new WorkerDownCommand(ID<worker_id_t>(worker->worker_id()));
    boost::mutex::scoped_lock lock(command_queue_mutex_);
    received_commands_.push_back(command);
    return;
  }

  total_bytes_received_ += bytes_transferred;

  size_t real_length = bytes_transferred + worker->connection()->existing_bytes();
  size_t len = EnqueueCommands(worker->connection()->read_buffer(), real_length);
  assert(real_length >= len);
  size_t remaining = (real_length - len);
  dbg(DBG_NET,
      "Scheduler received %i bytes from worker %i, enqueued %i bytes as commands, %i remaining.\n",
      bytes_transferred, worker->worker_id(), len, remaining);


  // This is the case when the buffer is full with incomplete message.
  // Then it cannot make progress anymore.
  if (remaining == worker->connection()->max_read_buffer_length()) {
    std::string err_msg = "ERROR: ";
    err_msg += "scheduler server buffer is full with incomplete command. ";
    err_msg += "You need to increase the SERVER_CON_BUF_SIZE in scheduler_server_connection.cc.\n";
    dbg(DBG_ERROR, "%s\n.", err_msg.c_str());
    exit(-1);
  }

  // This is for the case when the string buffer had an incomplete
  // command at its end. Copy the fragement of the command to the beginning
  // of the buffer, mark how many bytes are valid with existing_bytes.
  if (remaining > 0) {
    char* buffer = worker->connection()->read_buffer();
    memmove(buffer, (buffer + len), remaining);
    char* read_start_ptr = buffer + remaining;
    int read_buffer_length = worker->connection()->max_read_buffer_length() - remaining;

    worker->connection()->set_existing_bytes(remaining);
    boost::asio::async_read(*(worker->connection()->socket()),
                            boost::asio::buffer(read_start_ptr,
                                                read_buffer_length),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerServer::HandleRead,
                                        this,
                                        worker,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  } else {
    worker->connection()->set_existing_bytes(0);
    boost::asio::async_read(*(worker->connection()->socket()),
                            boost::asio::buffer(worker->connection()->read_buffer(),
                                                worker->connection()->max_read_buffer_length()),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerServer::HandleRead,
                                        this,
                                        worker,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  }
}


void SchedulerServer::HandleWrite(SchedulerWorker* worker,
                                  const boost::system::error_code& error,
                                  size_t bytes_transferred) {}


void SchedulerServer::Run() {
  dbg(DBG_NET, "Running the scheduler networking server.\n");
  Initialize();
  ListenForNewConnections();
  io_service_->run();
}

SchedulerWorkerList* SchedulerServer::workers() {
  return &workers_;
}

bool SchedulerServer::GetSchedulerWorkerById(SchedulerWorker*& worker, worker_id_t w_id) {
  SchedulerWorkerList::iterator iter = workers_.begin();
  for (; iter != workers_.end(); ++iter) {
    if ((*iter)->worker_id() == w_id) {
      worker = *iter;
      return true;
    }
  }
  return false;
}


bool SchedulerServer::RemoveWorker(worker_id_t worker_id) {
  SchedulerWorkerList::iterator iter = workers_.begin();
  for (; iter != workers_.end(); ++iter) {
    if ((*iter)->worker_id() == worker_id) {
      SchedulerWorker *worker = *iter;
      worker->MarkDead();
      workers_.erase(iter);
      delete worker;
      return true;
    }
  }

  dbg(DBG_WARN, "WARNING: could not find worker %lu to remove!\n", worker_id);
  exit(-1);
  return false;
}



size_t SchedulerServer::worker_num() {
  return workers_.size();
}

uint64_t SchedulerServer::total_bytes_sent() {
  return total_bytes_sent_;
}

uint64_t SchedulerServer::total_bytes_received() {
  return total_bytes_received_;
}

void SchedulerServer::set_bouncer_thread_active(bool flag) {
  bouncer_thread_active_ = flag;
}

void
SchedulerServer::set_worker_command_table(SchedulerCommand::PrototypeTable* cmt) {
  worker_command_table_ = cmt;
}

SchedulerWorker* SchedulerServer::AddWorker(SchedulerServerConnection* connection) { //NOLINT
  boost::mutex::scoped_lock lock(worker_mutex_);
  static worker_id_t workerIdentifier = NIMBUS_SCHEDULER_ID;
  workerIdentifier++;
  SchedulerWorker* worker = new SchedulerWorker(workerIdentifier,
                                                connection, NULL);
  workers_.push_back(worker);
  return worker;
}

void SchedulerServer::MarkWorkerDead(SchedulerWorker* worker) {
  worker->MarkDead();
}

}  // namespace nimbus
