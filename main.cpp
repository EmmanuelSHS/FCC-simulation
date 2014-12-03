#include "ThreeD.h"

int main(){
  ThreeD *md = new ThreeD();

  md -> init();
  md -> MD_loop();

  delete md;

return 0;
}
