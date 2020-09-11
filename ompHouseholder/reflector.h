
void GetReflector(int sz, const double * vector, double * reflector);
void ApplyReflector(int sz, const double * reflector, double * vec);
void ApplyReflectorArray(int sz, int block_sz, const double * reflector, int ldr, double * block, int ldblock);
