#ifndef __DATA_IMPL__
#define __DATA_IMPL__

#include "shared/nimbus.h"
#include "worker/data.h"

// A dirty version, for experiment.
// Data structure to describe partial norm.
class PartialNorm : public Data {
 public:
  explicit PartialNorm(double norm) {norm_ = norm;}
  virtual ~PartialNorm() {}
  virtual void Create() {printf("[QUH] Partial Norm created.\n");}
  virtual void Destroy() {}
  virtual Data* Clone() {
    printf("[QUH] Colne data.\n");
    return new PartialNorm(0);
  }
  virtual void Copy(Data* from) {
    PartialNorm* d = reinterpret_cast<PartialNorm*>(from);
    norm_ = d->norm_;
  }
  virtual bool Serialize(SerializedData* ser_data) {
    printf("[QUH] Serialization starts.\n");
    // TODO do serialization.
    double* buffer = new double;
    *buffer = norm_;
    ser_data->set_data_ptr((char*)buffer);
    ser_data->set_size(sizeof(*buffer));
    printf("[QUH] Serialization end%f.\n", norm_);
    return true;
  }
  virtual bool DeSerialize(const SerializedData& ser_data, Data** result) {
    printf("[QUH] DeSerialization starts.\n");
    double temp =  *(reinterpret_cast<double*>(ser_data.data_ptr_raw()));
    PartialNorm* partial_norm = new PartialNorm(temp);
    *result = partial_norm;
    printf("[QUH] DeSerialization end%f.\n", temp);
    return true;
  }
  // TODO moved to private member.
  double norm_;
};
#endif
