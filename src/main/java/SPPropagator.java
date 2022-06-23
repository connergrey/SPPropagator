import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.forces.ForceModel;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.DragSensitive;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.SolidTides;
import org.orekit.forces.gravity.ThirdBodyAttraction;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.forces.gravity.potential.TideSystem;
import org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient;
import org.orekit.forces.radiation.RadiationSensitive;
import org.orekit.forces.radiation.SolarRadiationPressure;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Transform;
import org.orekit.models.earth.atmosphere.DTM2000;
import org.orekit.models.earth.atmosphere.DTM2000InputParameters;
import org.orekit.models.earth.atmosphere.data.MarshallSolarActivityFutureEstimation;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.OrbitType;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.time.UT1Scale;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

import java.io.FileNotFoundException;

public class SPPropagator {

    public static void main(String[] arg) throws FileNotFoundException {

        // make sure you update orekit-data

        DataLoader loader = new DataLoader();
        loader.load(); // loads the orekit-data file to get constant parameters like tai-utc, etc.

        String solarInd = "mar2022f10_prd.txt"; // EXAMPLE... Please update to the most revelant solar weather file
        String spVectorFileName = "/Users/connergrey/Documents/SP VECTORS/vectors_22072/scratch/SP_VEC/33/33064"; // EXAMPLE... SP VECTOR FILE NAME

        // set constants
        double earthMu = Constants.EGM96_EARTH_MU;
        double erad = Constants.EGM96_EARTH_EQUATORIAL_RADIUS;
        TimeScale utc = TimeScalesFactory.getUTC();
        Frame eme2000 = FramesFactory.getEME2000();

        // read sp vector
        SPVecReader spvec = new SPVecReader();
        spvec.readSPVec(spVectorFileName);
        int degOrd = spvec.getDegOrd();
        double cR = spvec.getcR();
        double cD = spvec.getcD();

        // set initial state
        AbsoluteDate initialDate = spvec.getDate();
        TimeStampedPVCoordinates initial = spvec.getInitialState(eme2000);
        PVCoordinates pvCoordinatesEME2000 = new PVCoordinates(initial.getPosition(), initial.getVelocity());

        // Create orbit from initial condition
        CartesianOrbit orbit = new CartesianOrbit(pvCoordinatesEME2000, eme2000, initialDate, earthMu);

        // Create the integrator and propagator with desired models
        PropCreator propCreator = new PropCreator();
        NumericalPropagator propagator = PropCreator.createProp(orbit, earthMu, erad, utc, degOrd, cR, cD, solarInd);

        //propagate
        //double propTime = 3600;
        AbsoluteDate finalDate = new AbsoluteDate(2022,3,14,5,43,26.191,utc);
        SpacecraftState s = propagator.propagate(initialDate,finalDate); //output in meters in EME2000
        Frame teme = FramesFactory.getTEME();
        System.out.println(s.getPVCoordinates(teme).getPosition().scalarMultiply(1.0/1000));

    }

}
