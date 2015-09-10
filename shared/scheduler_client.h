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
  * Author: Omid Mashayekhi <omidm@stanford.edu>
  * Author: Philip Levis <pal@cs.stanford.edu>
  */

#ifndef NIMBUS_SHARED_SCHEDULER_CLIENT_H_
#define NIMBUS_SHARED_SCHEDULER_CLIENT_H_

#include <boost/thread.hpp>
#include <boost/asio.hpp>
#include <string>
#include <map>
#include <list>
#include <vector>
#include "shared/nimbus_types.h"
#include "shared/scheduler_command_include.h"
#include "shared/command_template.h"

namespace nimbus {

typedef uint32_t ConnectionId;

using boost::asio::ip::tcp;

/** Networking services for a Nimbus scheduler to talk to workers.
 *  All receive and send calls are non-blocking. Currently
 *  maintains a list of workers asynchronously: future versions
 *  will provide callbacks to notify scheduler of updates to the
 *  worker list.
 */
class SchedulerClient {
 public:
  SchedulerClient(std::string scheduler_ip, port_t scheduler_port);
  virtual ~SchedulerClient();

  /** Start the client, does not return. A running client starts a connection
   * to server, spools received messages into its message queue, and services
   * send requests. */
  virtual void Run();

  /** Pull incoming commands from the received queue. Puts at most
   *  maxCommands into storage, returning true if it placed one or more.
   *  Returns false if no commands were placed in storage. */
  virtual bool ReceiveCommands(SchedulerCommandList* storage,
                               size_t maxCommands);

  /** Send command to server. Returns immediately and processes the send
   * asynchronously.*/
  virtual void SendCommand(SchedulerCommand* command);

  /** Send commands to server. Returns immediately and processes the send
   * asynchronously. */
  virtual void SendCommands(SchedulerCommandList* commands);

  void set_scheduler_command_table(SchedulerCommand::PrototypeTable* cmt);

  void set_command_processor_mutex(boost::recursive_mutex *mutex);

  void set_command_processor_cond(boost::condition_variable_any *cond);

  void PushCommandToTheQueue(SchedulerCommand *command);

 private:
  std::string scheduler_ip_;
  port_t scheduler_port_;
  boost::asio::io_service* io_service_;
  tcp::socket* socket_;
  SchedulerCommand::PrototypeTable* scheduler_command_table_;
  char* read_buffer_;
  uint32_t max_read_buffer_length_;
  uint32_t existing_bytes_;
  SchedulerCommandList received_commands_;
  boost::mutex send_command_mutex_;

  boost::recursive_mutex *command_processor_mutex_;
  boost::condition_variable_any *command_processor_cond_;

  boost::thread* command_template_thread_;
  typedef std::map<std::string, CommandTemplate*> TemplateMap;
  TemplateMap template_map_;
  bool filling_template_;
  std::string template_name_in_progress_;

  class CommandTemplateSeed {
    public:
      CommandTemplateSeed(CommandTemplate *command_template,
                          const std::vector<job_id_t>& inner_job_ids,
                          const std::vector<job_id_t>& outer_job_ids,
                          const std::vector<Parameter>& parameters,
                          const std::vector<physical_data_id_t>& physical_ids)
        : command_template_(command_template),
          inner_job_ids_(inner_job_ids),
          outer_job_ids_(outer_job_ids),
          parameters_(parameters),
          physical_ids_(physical_ids) {}
      ~CommandTemplateSeed() {}

      CommandTemplate *command_template_;
      std::vector<job_id_t> inner_job_ids_;
      std::vector<job_id_t> outer_job_ids_;
      std::vector<Parameter> parameters_;
      std::vector<physical_data_id_t> physical_ids_;
  };

  typedef std::list<CommandTemplateSeed*> CommandTemplateSeedList;
  CommandTemplateSeedList command_template_seeds_;
  boost::mutex command_template_seeds_mutex_;
  boost::condition_variable_any command_template_seeds_cond_;

  virtual void CommandTemplateThread();

  /** Create client socket, set up networking and state. */
  virtual bool Initialize();

  /** Create new connection to the server. */
  void CreateNewConnections();

  /** Asynchronous callback to read data. */
  virtual void HandleRead(const boost::system::error_code& error,
                          size_t bytes_transferred);

  /** Asynchronous calback for writing data. */
  virtual void HandleWrite(const boost::system::error_code& error,
                           size_t bytes_transferred);

  /** Call to take a buffer from the network of size, parse them
      into commands and put them on the pending command queue.
      Return value is how many bytes were read/parsed, this can
      be less than size (for example, if the buffer end has an
      incomplete command. */
  virtual size_t EnqueueCommands(char* buffer, size_t size);

  virtual void HandleStartCommandTemplateCommand(StartCommandTemplateCommand *cm);

  virtual void HandleEndCommandTemplateCommand(EndCommandTemplateCommand *cm);

  virtual void HandleSpawnCommandTemplateCommand(SpawnCommandTemplateCommand *cm);
};

}  // namespace nimbus

#endif  // NIMBUS_SHARED_SCHEDULER_CLIENT_H_
