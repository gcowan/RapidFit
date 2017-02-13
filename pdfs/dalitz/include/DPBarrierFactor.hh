#ifndef DP_BARRIER_FUNCTION
#define DP_BARRIER_FUNCTION
#include <memory>
class DPBarrierFactor
{
	public:
		DPBarrierFactor();
		DPBarrierFactor(const unsigned spin, const double radius, const double p0);
		DPBarrierFactor(const DPBarrierFactor& other) : spin(other.spin), function(other.function), precalcFF(other.precalcFF), radius(other.radius) {}
		~DPBarrierFactor() {}
		double barrier(const double p) const;
		double barrier_sq(const double p) const;
		int getspin() const {return spin;}
		void setparameters(const double _radius, const double p0);
		static double FunctionL0(const double z);
		static double FunctionL1(const double z);
		static double FunctionL2(const double z);
		static double FunctionL3(const double z);
	private:
		double (*function)(const double z);
		double radius;  // Blatt-Weiskopf radius
		double precalcFF; // Pre-calculated form factor using p0
		int spin;
};

#endif
