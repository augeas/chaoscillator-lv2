#include <math.h>
#include <stdlib.h>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <lv2/lv2plug.in/ns/lv2core/lv2.h>

#define LORENZ_URI "https://augeas.github.io/chaoscillator-lv2/lorenz"

using namespace boost::numeric::odeint;

enum PortIndex {
    LORENZ_SIGMA,
    LORENZ_RHO,
    LORENZ_BETA,
    LORENZ_X,
    LORENZ_Y,
    LORENZ_Z
};

typedef struct {
    const float* sigma;
    const float* rho;
    const float* beta;
    float* x;
    float* y;
    float* z;
    float xx;
    float yy;
    float zz;
} Lorenz;

static LV2_Handle
instantiate(const LV2_Descriptor* descriptor, double rate, const char* bundle_path,
    const LV2_Feature* const* features) {
        Lorenz* lorenz = (Lorenz*)malloc(sizeof(Lorenz));
        return (LV2_Handle)lorenz;
}

static void
connect_port(LV2_Handle instance, uint32_t port, void* data) {
    Lorenz* lorenz = (Lorenz*)instance;
    PortIndex pindex = (PortIndex) port;
    
    switch (pindex) {
        case LORENZ_SIGMA:
            lorenz->sigma = (const float*)data;
            break;
        case LORENZ_RHO:
            lorenz->rho = (const float*)data;
            break;
        case LORENZ_BETA:
            lorenz->beta = (float*)data;
            break;
        case LORENZ_X:
            lorenz->x = (float*)data;
            break;
        case LORENZ_Y:
            lorenz->y = (float*)data;
            break;
        case LORENZ_Z:
            lorenz->z = (float*)data;
            break;  
    }
}

static void
activate(LV2_Handle instance)
{
    Lorenz* lorenz = (Lorenz*)instance;
    lorenz->xx = 0.1;
    lorenz->yy = 0.1;
    lorenz->zz = 0.1;
}

typedef std::vector<float> state_type;

static void
run(LV2_Handle instance, uint32_t n_samples)
{
        state_type vec(3);
    
        Lorenz* lorenz = (Lorenz*)instance;

        const float* const sigma = lorenz->sigma;
        const float* const rho = lorenz->rho;
        const float* const beta = lorenz->beta;
        
        float* const x = lorenz->x;
        float* const y = lorenz->y;
        float* const z = lorenz->z;
        
        vec[0] = lorenz->xx;
        vec[1] = lorenz->yy;
        vec[2] = lorenz->zz;
        
        runge_kutta4<state_type> rk4;
        
        for (uint32_t pos = 0; pos < n_samples; pos++) {
            rk4.do_step([&sigma, &rho, &beta, &pos] (state_type &v , state_type &dvdt , double t) {
                dvdt[0] = sigma[pos] * (v[1] - v[0]);
                dvdt[1] = rho[pos] * v[0] - v[1] - v[0] * v[2];
                dvdt[2] = v[0]*v[1] - beta[pos] * v[2];
            }, vec ,0.0 ,0.01);
            x[pos] = vec[0];
            y[pos] = vec[1];
            z[pos] = vec[2];
        }
        
        lorenz->xx = vec[0];
        lorenz->yy = vec[1];
        lorenz->zz = vec[2];
}

static void
deactivate(LV2_Handle instance)
{
}

static void
cleanup(LV2_Handle instance)
{
    free(instance);
}

static const void*
extension_data(const char* uri)
{
    return NULL;
}

static const LV2_Descriptor descriptor = {
    LORENZ_URI,
    instantiate,
    connect_port,
    activate,
    run,
    deactivate,
    cleanup,
    extension_data
};

LV2_SYMBOL_EXPORT
const LV2_Descriptor*
lv2_descriptor(uint32_t index)
{
    switch (index) {
        case 0:  return &descriptor;
        default: return NULL;
    }
}
