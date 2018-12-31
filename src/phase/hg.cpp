/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/medium.h>
#include <mitsuba/render/phase.h>
#include <mitsuba/render/volume.h>
#include <mitsuba/render/sampler.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/frame.h>

MTS_NAMESPACE_BEGIN

/*!\plugin{hg}{Henyey-Greenstein phase function}
 * \order{2}
 * \parameters{
 *     \parameter{g}{\Float}{
 *       This parameter must be somewhere in the range $-1$ to $1$
 *       (but not equal to $-1$ or $1$). It denotes the \emph{mean cosine}
 *       of scattering interactions. A value greater than zero indicates that
 *       medium interactions predominantly scatter incident light into a similar
 *       direction (i.e. the medium is \emph{forward-scattering}), whereas
 *       values smaller than zero cause the medium to be
 *       scatter more light in the opposite direction.
 *     }
 * }
 * This plugin implements the phase function model proposed by
 * Henyey and Greenstein \cite{Henyey1941Diffuse}. It is
 * parameterizable from backward- ($g<0$) through
 * isotropic- ($g=0$) to forward ($g>0$) scattering.
 */
class HGPhaseFunction : public PhaseFunction {
public:

	HGPhaseFunction(const Properties &props)
		: PhaseFunction(props) {
	}
	
	/* Unserialize from a binary data stream */
	HGPhaseFunction(Stream *stream, InstanceManager *manager)
		: PhaseFunction(stream, manager) {
		m_g = static_cast<VolumeDataSource *>(manager->getInstance(stream));
		configure();
	}

	virtual ~HGPhaseFunction() { }

	void serialize(Stream *stream, InstanceManager *manager) const {
		PhaseFunction::serialize(stream, manager);
		manager->serialize(stream, m_g.get());
	}

	void configure() {
		PhaseFunction::configure();
		if (m_g.get() == NULL)
			Log(EError, "No g specified!");
		m_type = EAngleDependence;
	}

	void addChild(const std::string &name, ConfigurableObject *child) {
		if (child->getClass()->derivesFrom(MTS_CLASS(VolumeDataSource))) {
			Assert(name == "g");
			m_g = static_cast<VolumeDataSource *>(child);
		}
	}
	
	
	inline Float sample(PhaseFunctionSamplingRecord &pRec,
			Sampler *sampler) const {
		Point2 sample(sampler->next2D());
		Point p = pRec.mRec.p;
		Float g = m_g->lookupFloat(p);
		Float cosTheta;
		if (std::abs(g) < Epsilon) {
			cosTheta = 1 - 2*sample.x;
		} else {
			Float sqrTerm = (1 - g * g) / (1 - g + 2 * g * sample.x);
			cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
		}

		Float sinTheta = math::safe_sqrt(1.0f-cosTheta*cosTheta),
			  sinPhi, cosPhi;

		math::sincos(2*M_PI*sample.y, &sinPhi, &cosPhi);

		pRec.wo = Frame(-pRec.wi).toWorld(Vector(
			sinTheta * cosPhi,
			sinTheta * sinPhi,
			cosTheta
		));
		return 1.0f;
	}

	Float sample(PhaseFunctionSamplingRecord &pRec,
			Float &pdf, Sampler *sampler) const {
		HGPhaseFunction::sample(pRec, sampler);
		pdf = HGPhaseFunction::eval(pRec);
		return 1.0f;
	}

	Float eval(const PhaseFunctionSamplingRecord &pRec) const {
		Point p = pRec.mRec.p;
		Float g = m_g->lookupFloat(p);
		Float temp = 1.0f + g*g + 2.0f * g * dot(pRec.wi, pRec.wo);
		return INV_FOURPI * (1 - g*g) / (temp * std::sqrt(temp));
	}

	Float getMeanCosine() const {
		return 0;
	}
	
	Float getMeanCosine(Point p) const {
		return m_g->lookupFloat(p);
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "HGPhaseFunction[g=" << m_g.toString() << "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
	ref<VolumeDataSource> m_g;
};

MTS_IMPLEMENT_CLASS_S(HGPhaseFunction, false, PhaseFunction)
MTS_EXPORT_PLUGIN(HGPhaseFunction, "Henyey-Greenstein phase function");
MTS_NAMESPACE_END
