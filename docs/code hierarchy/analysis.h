
class analysis {
public:
    analysis(circuit * c);
    enum TYPE_t {
        DC      = 0,
        TRAN_FE = 1,
        TRAN_BE = 2,
        TRAN_TR = 3
    };
    virtual void plotnv(matlab * m, int node_name) = 0;
    virtual void printnv(int node_name) = 0;
private:
    static constexpr double precision = 1.0e-9;
    circuit * c;
    std::vector<circuit::linelem*> itrelems;
};
