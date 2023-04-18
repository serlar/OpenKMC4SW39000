extern "C"{
    void init_prof(const char *path);
    void start_prof();
    void pause_prof();
    void init_prof_rank(int iid);
    void stop_prof();
}

