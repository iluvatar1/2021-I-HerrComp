#include <iostream>
#include <string>

int main(int argc, char **argv)
{
  std::string firstname, lastname, fullname;
  firstname = "William";
  lastname = "Oquendo";
  fullname = firstname + " " + lastname;

  std::cout << "Hola: " << fullname << "\n";
  std::cout << fullname[5] << "\n";
  std::cout << fullname.size() << "\n";
  std::cout << fullname.max_size() << "\n";

  std::string a = "7007.0909979";
  double x = std::stod(a);
  std::cout << x*2 << std::endl;

  return 0;
}

