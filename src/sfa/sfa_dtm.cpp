#include "sfa.h"
#include "utility/logger.h"
#include <algorithm>

namespace SFA {
    cvec2 SFA::DTM(const dvec2& p) {
        static const double dp = ps[1] - ps[0];
        static const double invdp = 1. / dp;

        double l = length(p);
        double di = (l - ps[0]) * invdp;
        int i = int(di);
        i = std::max(0, std::min(i, int(ps.size()-2)));
        int j = i + 1;
        LOG_DEBUG(std::to_string(ps.size()));
        LOG_DEBUG(std::to_string(j));
        if (di == double(i))
            return dtm_el[i] * p;
            


        double f = (l - ps[i]) * invdp;
        complex dtm = dtm_el[j] * f + (1. - f) * dtm_el[i];

        return dtm * p;
    }
}