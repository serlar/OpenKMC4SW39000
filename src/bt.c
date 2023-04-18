#include <execinfo.h>
#include <stdio.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <ucontext.h>
static int out_fd;
static void **out_buffer;
static int buffer_ptr;
static int buffer_cnt = 1024 * 1024 * 16;
static timer_t timerid;
static jmp_buf portal;
static struct sigaction act;
static int in_alrm = 0;
// #ifdef __OPEN64__
// #include <swlu.h>
// //int swlu_backtrace(void **buffer, int size);
// #define backtrace(x, y) swlu_backtrace_context(context, x, y)
// #endif
void flush_buffer(){
  write(out_fd, out_buffer, buffer_ptr * sizeof(void*));
  buffer_ptr = 0;
}
void sig_segv(int a){
  longjmp(portal, 1);
}
void sig_alrm_prof(int a, siginfo_t *info, ucontext_t *context){
  void **buffer = out_buffer + buffer_ptr;
  //#ifndef __OPEN64__
  if (in_alrm) return;
  in_alrm = 1;
  signal(SIGSEGV, sig_segv);
  //#endif
  int error = setjmp(portal);
  int ntrace;
  if (!error){
    ntrace = backtrace(buffer, 100);
  } else {
    ntrace = 0;
  }
  buffer[ntrace] = 0;
  buffer_ptr += ntrace + 1;
  if (buffer_ptr + 101 >= buffer_cnt) {
    flush_buffer();
  }
  signal(SIGSEGV, SIG_DFL);
  in_alrm = 0;
}
void sig_alrm_stub(int a, siginfo_t *info, ucontext_t *context){
}

void init_prof(const char *path){
  out_buffer = malloc(buffer_cnt * sizeof(void *));
  buffer_ptr = 0;
  out_fd = open(path, O_WRONLY | O_CREAT | O_TRUNC, 0755);
  sigaction(SIGALRM, NULL, &act);
  //act.sa_flags = SA_ONSTACK | SA_NODEFER;
  act.sa_flags |= SA_NODEFER | SA_ONSTACK | SA_SIGINFO;
  //act.sa_flags &= ~SA_NODEFER;
  act.sa_sigaction = (void*)sig_alrm_stub;
  sigaction(SIGALRM, &act, NULL);

  struct sigevent sev;
  struct itimerspec its;

  timer_create(CLOCK_REALTIME, NULL, &timerid);
  its.it_value.tv_sec = 0;
  its.it_value.tv_nsec = 1000000;
  its.it_interval.tv_sec = 0;
  its.it_interval.tv_nsec = 1000000;
  timer_settime(timerid, 0, &its, NULL);

  stack_t sigstk;
  sigstk.ss_size = 0;
  sigstk.ss_flags = 0;
  sigstk.ss_sp = malloc (SIGSTKSZ);
  if (sigstk.ss_sp != NULL) {
    sigstk.ss_size = SIGSTKSZ;
    if (sigaltstack (&sigstk, 0) < 0) {
      sigstk.ss_size = 0;
      free (sigstk.ss_sp);
      perror("sigaltstack error");
    }
  } else {
    printf ("malloc (SIGSTKSZ) failed!\n");
  }

  //signal(SIGALRM, sig_alrm_stub);
}
void start_prof(){
  //signal(SIGALRM, sig_alrm_prof);
  act.sa_sigaction = (void*)sig_alrm_prof;
  sigaction(SIGALRM, &act, NULL);

}
void pause_prof(){
  //signal(SIGALRM, sig_alrm_stub);
  act.sa_sigaction = (void*)sig_alrm_stub;
  sigaction(SIGALRM, &act, NULL);
}
void init_prof_rank(int iid) {
  char path[100];
  strcpy(path, "backtrace.");
  int ptr = strlen(path);
  path[ptr + 0] = '0' + iid / 100000;
  path[ptr + 1] = '0' + iid % 100000 / 10000;
  path[ptr + 2] = '0' + iid % 10000 / 1000;
  path[ptr + 3] = '0' + iid % 1000 / 100;
  path[ptr + 4] = '0' + iid % 100 / 10;
  path[ptr + 5] = '0' + iid % 10;
  path[ptr + 6] = 0;
  strcat(path, ".bin");
  init_prof(path);
}
void stop_prof(){
  timer_delete(timerid);
  flush_buffer();
}
void init_prof_(int *id) {init_prof_rank(*id);}
void stop_prof_() __attribute__((alias("stop_prof")));
void start_prof_() __attribute__((alias("start_prof")));
void pause_prof_() __attribute__((alias("pause_prof")));

#ifdef TEST
int main(){
  int id = 0;
  init_prof_(&id);
  start_prof();
  for (int i = 0; i < 100000; i ++){
    printf("%d\n", i);
  }
  pause_prof();
  stop_prof();
  return 0;
}
#endif
