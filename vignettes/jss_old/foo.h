// test document embedding

// basic setter/getter for scalar double
class foo {
 private:
  double x;
 public:
  foo(double _x);
  void set_x(double _x);
  double get_x();
};

// constructor
inline foo::foo(double _x) {
  set_x(_x);
}

// getter
double get_x() {
  return x;
}

// setter
void set_x(double _x) {
  x = _x;
  return;
}
