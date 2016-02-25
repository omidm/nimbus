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
  * Client (worker) side interface of the Nimbus scheduler protocol.
  *
  * This interface is a simple message passing interface to the scheduler.
  * It's responsible for transforming commands between object and wire
  * formats. On the reception side, this means instantiating the proper
  * command subclass using a factory.

  * Calls to the network should be non-blocking. In the case of
  * sending, this means we need to have error return values for
  * whether the send was successful. It should also properly handle no
  * messages to receive through non-blocking read calls.
  *
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */


#include <boost/tokenizer.hpp>
#include <string>
#include <sstream>
#include <iostream>  // NOLINT
#include "src/shared/dbg.h"
#include "src/shared/scheduler_client.h"


using boost::asio::ip::tcp;

#define CLIENT_TCP_SEND_BUF_SIZE 134217728  // 128MB
#define CLIENT_TCP_RECEIVE_BUF_SIZE 134217728  // 128MB
#define CLIENT_BUF_SIZE 40960000
#define SEPARATE_COMMAND_TEMPLATE_THREAD true

namespace nimbus {


SchedulerClient::SchedulerClient(std::string scheduler_ip,
                                 port_t scheduler_port)
  : scheduler_ip_(scheduler_ip),
    scheduler_port_(scheduler_port) {
      filling_template_ = false;
      execution_template_active_ = true;
}

SchedulerClient::~SchedulerClient() {
  delete socket_;
  delete read_buffer_;
  delete io_service_;
}

bool SchedulerClient::Initialize() {
  io_service_ = new boost::asio::io_service();
  socket_ = new tcp::socket(*io_service_);
  read_buffer_ = new char[CLIENT_BUF_SIZE];;
  max_read_buffer_length_ = CLIENT_BUF_SIZE;
  existing_bytes_ = 0;
  return true;
}

bool SchedulerClient::ReceiveCommands(SchedulerCommandList* storage,
                                      size_t maxCommands) {
  storage->clear();
  boost::unique_lock<boost::recursive_mutex> lock(*command_processor_mutex_);
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

void SchedulerClient::SendCommand(SchedulerCommand* command) {
  std::string data = command->ToNetworkData();
  SchedulerCommand::length_field_t len;
  len = htonl((uint32_t)data.length() + sizeof(len));
  std::string msg;
  msg.append((const char*)&len, sizeof(len));
  msg.append(data.c_str(), data.length());

  boost::mutex::scoped_lock lock(send_command_mutex_);
  // dbg(DBG_NET, "Sending command %s.\n", command->ToString().c_str());
  boost::system::error_code ignored_error;
  // Why are we IGNORING ERRORS!??!?!?
  boost::asio::write(*socket_, boost::asio::buffer(msg),
                     boost::asio::transfer_all(), ignored_error);
}

void SchedulerClient::SendCommands(SchedulerCommandList* commands) {
  SchedulerCommandList::iterator iter = commands->begin();
  for (; iter != commands->end(); ++iter) {
    SchedulerCommand* command = *iter;
    SendCommand(command);
  }
}

using boost::tokenizer;
using boost::char_separator;

  /** Reads commands from buffer with size bytes. Puts commands into
   * internal command list for later reads. Returns the number of
   * bytes read into commands. If this is less than size, the
   * caller must save the unread bytes. */
size_t SchedulerClient::EnqueueCommands(char* buffer, size_t size) {
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
                                                          scheduler_command_table_,
                                                          command)) {
        if (execution_template_active_) {
          dbg(DBG_NET, "Enqueueing command %s.\n", command->ToString().c_str());
          boost::unique_lock<boost::recursive_mutex> lock(*command_processor_mutex_);
          received_commands_.push_back(command);
          command_processor_cond_->notify_all();
        } else {
          switch (command->type()) {
            case SchedulerCommand::START_COMMAND_TEMPLATE:
              HandleStartCommandTemplateCommand(
                  reinterpret_cast<StartCommandTemplateCommand*>(command));
              break;
            case SchedulerCommand::END_COMMAND_TEMPLATE:
              HandleEndCommandTemplateCommand(
                  reinterpret_cast<EndCommandTemplateCommand*>(command));
              break;
            case SchedulerCommand::SPAWN_COMMAND_TEMPLATE:
              HandleSpawnCommandTemplateCommand(
                  reinterpret_cast<SpawnCommandTemplateCommand*>(command));
              break;
            default:
              if (filling_template_) {
                TemplateMap::iterator iter = template_map_.find(template_name_in_progress_);
                assert(iter != template_map_.end());
                CommandTemplate *ct = iter->second;
                switch (command->type()) {
                  case SchedulerCommand::EXECUTE_COMPUTE:
                    ct->AddComputeJobCommand(reinterpret_cast<ComputeJobCommand*>(command));
                    break;
                  case SchedulerCommand::LOCAL_COPY:
                    ct->AddLocalCopyCommand(reinterpret_cast<LocalCopyCommand*>(command));
                    break;
                  case SchedulerCommand::REMOTE_SEND:
                    ct->AddRemoteCopySendCommand(reinterpret_cast<RemoteCopySendCommand*>(command));
                    break;
                  case SchedulerCommand::REMOTE_RECEIVE:
                    ct->AddRemoteCopyReceiveCommand(reinterpret_cast<RemoteCopyReceiveCommand*>(command)); // NOLINT
                    break;
                  case SchedulerCommand::MEGA_RCR:
                    ct->AddMegaRCRCommand(reinterpret_cast<MegaRCRCommand*>(command)); // NOLINT
                    break;
                  default:
                    assert(false);
                    break;
                }
              }
              dbg(DBG_NET, "Enqueueing command %s.\n", command->ToString().c_str());
              boost::unique_lock<boost::recursive_mutex> lock(*command_processor_mutex_);
              received_commands_.push_back(command);
              command_processor_cond_->notify_all();
              break;
          }
        }
      } else {
        dbg(DBG_NET, "Ignored unknown command: %s.\n", input.c_str());
        assert(false);
      }
      dbg(DBG_NET, "Read %i bytes of %i available.\n", len, bytes_remaining);
      bytes_remaining -= len;
      total_read += len;
      read_ptr += len;
    }
  }

  return total_read;
}

void SchedulerClient::HandleStartCommandTemplateCommand(StartCommandTemplateCommand *cm) {
  std::string key = cm->command_template_name();
  assert(template_map_.find(key) == template_map_.end());
  CommandTemplate *ct = new CommandTemplate(key,
                                            cm->inner_job_ids(),
                                            cm->outer_job_ids(),
                                            cm->phy_ids());
  template_map_[key] = ct;
  filling_template_ = true;
  template_name_in_progress_ = key;
}

void SchedulerClient::HandleEndCommandTemplateCommand(EndCommandTemplateCommand *cm) {
  assert(filling_template_);
  std::string key = cm->command_template_name();
  assert(template_name_in_progress_ == key);
  TemplateMap::iterator iter = template_map_.find(key);
  assert(iter != template_map_.end());
  iter->second->Finalize();
  filling_template_ = false;
}

void SchedulerClient::HandleSpawnCommandTemplateCommand(SpawnCommandTemplateCommand *cm) {
  std::string key = cm->command_template_name();
  TemplateMap::iterator iter = template_map_.find(key);
  assert(iter != template_map_.end());
  if (SEPARATE_COMMAND_TEMPLATE_THREAD) {
    CommandTemplateSeed *cts =
      new CommandTemplateSeed(iter->second,
                              cm->inner_job_ids(),
                              cm->outer_job_ids(),
                              cm->extra_dependency(),
                              cm->parameters(),
                              cm->phy_ids());
    {
      boost::unique_lock<boost::mutex> command_template_seeds_lock(command_template_seeds_mutex_);
      command_template_seeds_.push_back(cts);
      command_template_seeds_cond_.notify_all();
    }
  } else {
  iter->second->Instantiate(cm->inner_job_ids(),
                            cm->outer_job_ids(),
                            cm->extra_dependency(),
                            cm->parameters(),
                            cm->phy_ids(),
                            this);
  }
}

void SchedulerClient::HandleRead(const boost::system::error_code& error,
                                 size_t bytes_transferred) {
  if (error) {
    dbg(DBG_NET|DBG_ERROR,
        "Error %s receiving %i bytes from scheduler.\n",
        error.message().c_str(), bytes_transferred);
    return;
  }

  size_t real_length = bytes_transferred + existing_bytes_;
  size_t len = EnqueueCommands(read_buffer_, real_length);
  assert(real_length >= len);
  size_t remaining = (real_length - len);
  dbg(DBG_NET,
      "Worker received %i bytes from scheduler, enqueued %i bytes as commands, %i remaining.\n",
      bytes_transferred, len, remaining);


  // This is the case when the buffer is full with incomplete message.
  // Then it cannot make progress anymore.
  if (remaining == max_read_buffer_length_) {
    std::string err_msg = "ERROR: ";
    err_msg += "scheduler client buffer is full with incomplete command. ";
    err_msg += "You need to increase the CLIENT_BUF_SIZE in scheduler_client.cc.\n";
    dbg(DBG_ERROR, "%s\n.", err_msg.c_str());
    exit(-1);
  }

  // This is for the case when the string buffer had an incomplete
  // command at its end. Copy the fragement of the command to the beginning
  // of the buffer, mark how many bytes are valid with existing_bytes.
  if (remaining > 0) {
    char* buffer = read_buffer_;
    memmove(buffer, (buffer + len), remaining);
    char* read_start_ptr = buffer + remaining;
    int read_buffer_length = max_read_buffer_length_ - remaining;

    existing_bytes_ = remaining;
    boost::asio::async_read(*socket_,
                            boost::asio::buffer(read_start_ptr,
                                                read_buffer_length),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerClient::HandleRead,
                                        this,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  } else {
    existing_bytes_ = 0;
    boost::asio::async_read(*socket_,
                            boost::asio::buffer(read_buffer_,
                                                max_read_buffer_length_),
                            boost::asio::transfer_at_least(1),
                            boost::bind(&SchedulerClient::HandleRead,
                                        this,
                                        boost::asio::placeholders::error,
                                        boost::asio::placeholders::bytes_transferred));
  }
}


void SchedulerClient::HandleWrite(const boost::system::error_code& error,
                                  size_t bytes_transferred) {}


void SchedulerClient::CreateNewConnections() {
  dbg(DBG_NET, "Attempting to open a connection to server.\n");
  tcp::resolver resolver(*io_service_);
  tcp::resolver::query query(scheduler_ip_,
      boost::to_string(scheduler_port_));
  tcp::resolver::iterator iterator = resolver.resolve(query);

  while (true) {
    boost::system::error_code error = boost::asio::error::host_not_found;
    socket_->connect(*iterator, error);
    if (error) {
      dbg(DBG_ERROR, "Could not connect to server: %s , retrying in 1 second ... \n", error.message().c_str()); // NOLINT
      delete socket_;
      socket_ = new tcp::socket(*io_service_);
      sleep(1);
    } else {
      dbg(DBG_NET, "Successfully connected to the server.\n");
      break;
    }
  }

  // Set the tcp send and receive buf size.
  // Note: you may have to increase the OS limits first.
  // Look at the nimbus/scripts/configure_tcp.sh for help.
  boost::asio::socket_base::send_buffer_size s_option(CLIENT_TCP_SEND_BUF_SIZE);
  boost::asio::socket_base::receive_buffer_size r_option(CLIENT_TCP_RECEIVE_BUF_SIZE);
  socket_->set_option(s_option);
  socket_->set_option(r_option);
  // Turn of Nagle algorithm.
  boost::asio::ip::tcp::no_delay nd_option(TCP_NODELAY_OPTION);
  socket_->set_option(nd_option);
  boost::asio::async_read(*socket_,
                          boost::asio::buffer(read_buffer_,
                                              max_read_buffer_length_),
                          boost::asio::transfer_at_least(1),
                          boost::bind(&SchedulerClient::HandleRead,
                                      this,
                                      boost::asio::placeholders::error,
                                      boost::asio::placeholders::bytes_transferred));
}


void SchedulerClient::Run() {
  dbg(DBG_NET, "Running the scheduler client.\n");
  Initialize();
  if (SEPARATE_COMMAND_TEMPLATE_THREAD) {
    command_template_thread_ =
      new boost::thread(boost::bind(&SchedulerClient::CommandTemplateThread, this));
  }
  CreateNewConnections();
  io_service_->run();
}

void SchedulerClient::CommandTemplateThread() {
  while (true) {
    CommandTemplateSeed *cts;
    {
      boost::unique_lock<boost::mutex> command_template_seeds_lock(command_template_seeds_mutex_);

      while (command_template_seeds_.size() == 0) {
        command_template_seeds_cond_.wait(command_template_seeds_lock);
      }

      CommandTemplateSeedList::iterator iter = command_template_seeds_.begin();
      cts = *iter;
      command_template_seeds_.erase(iter);
    }

    cts->command_template_->Instantiate(cts->inner_job_ids_,
                                        cts->outer_job_ids_,
                                        cts->extra_dependency_,
                                        cts->parameters_,
                                        cts->physical_ids_,
                                        this);
    delete cts;
  }
}

void
SchedulerClient::set_scheduler_command_table(SchedulerCommand::PrototypeTable* cmt) {
  scheduler_command_table_ = cmt;
}

void SchedulerClient::set_command_processor_mutex(boost::recursive_mutex *mutex) {
  command_processor_mutex_ = mutex;
}

void SchedulerClient::set_command_processor_cond(boost::condition_variable_any *cond) {
  command_processor_cond_ = cond;
}

void SchedulerClient::PushCommandToTheQueue(SchedulerCommand *command) {
  boost::unique_lock<boost::recursive_mutex> lock(*command_processor_mutex_);
  received_commands_.push_back(command);
  command_processor_cond_->notify_all();
}

void SchedulerClient::set_execution_template_active(bool flag) {
  execution_template_active_ = flag;
}


}  // namespace nimbus
