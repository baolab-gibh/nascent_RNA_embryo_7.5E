#include "SamAuxTag.h"

template <class T>
SamAuxTag<T>::SamAuxTag(
  const std::string& n, const std::string& t, const T& v
) : name(n), type(t), value(v)
{ }

template <class T>
SamAuxTag<T>::~SamAuxTag() = default;

template <class T>
const std::string& SamAuxTag<T>::SamAuxTag::get_name() const { return this->name; }

template <class T>
const std::string& SamAuxTag<T>::SamAuxTag::get_type() const { return this->type; }

template <class T>
const T& SamAuxTag<T>::SamAuxTag::get_value() const { return this->value; }

// Format the current object into a string.
template <>
const std::string SamAuxTag<std::string>::to_string() const {
  return this->name + ":" + this->type + ":" + this->value;
}

template <class T>
const std::string SamAuxTag<T>::to_string() const {
  return this->name + ":" + this->type + ":" + std::to_string(this->value);
}

// Cast the current object's value to the corresponding type.
// int.
template<>
int SamAuxTag<std::string>::get_val_int() const {
  if (this->value.empty()) {
    return 0;
  }

  if (this->type != "i") {
    std::cerr << "Error: Invalid type: " << this->type << std::endl;
    return 0;
  }

  return std::stoi(this->value);
}

// float
template<>
float SamAuxTag<std::string>::get_val_float() const {
  if (this->value.empty()) {
    return 0;
  }

  if (this->type != "f" || this->type != "d") {
    std::cerr << "Error: Invalid type: " << this->type << std::endl;
    return 0;
  }

  return std::stof(this->value);
}

// double
template<>
double SamAuxTag<std::string>::get_val_double() const { return SamAuxTag<std::string>::get_val_float(); }
