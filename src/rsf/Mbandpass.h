#ifndef H_BANDPASS
#define H_BANDPASS
void filter(float *trace, int nt, float dt, float flo, float fhi, bool phase, bool verb);
int main2(int argc, char* argv[]);
#endif 
