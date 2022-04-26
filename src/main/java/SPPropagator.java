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

        String solarInd = "feb2022f10_prd.txt"; // EXAMPLE... Please update to the most revelant solar weather file
        String spVectorFileName = "33063_2022056"; // EXAMPLE... SP VECTOR FILE NAME

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
        NumericalPropagator propagator = createProp(orbit, earthMu, erad, utc, degOrd, cR, cD, solarInd);

        //propagate
        double propTime = 3600;
        SpacecraftState s = propagator.propagate(initialDate,initialDate.shiftedBy(propTime)); //output in meters in EME2000

    }

    public static NumericalPropagator createProp(CartesianOrbit orbit, double earthMu, double erad, TimeScale utc,
                                                 int degOrd, double cR, double dragCoeff, String solarInd){


        IERSConventions iers = IERSConventions.IERS_2010;
        Frame itrf = FramesFactory.getITRF(iers,false);
        UT1Scale ut1 = TimeScalesFactory.getUT1(iers,false);
        OneAxisEllipsoid earth = new OneAxisEllipsoid(erad, Constants.WGS84_EARTH_FLATTENING, itrf);

        CelestialBody sun = CelestialBodyFactory.getSun();
        CelestialBody moon = CelestialBodyFactory.getMoon();

        double area = 1; // SRP and DRAG area considered equal
        int degree = degOrd;
        int order = degOrd;
        int mass = 1;

        SpacecraftState initialState = new SpacecraftState(orbit, mass);

        // Setup the integrator and propagator
        double dPos = 1e-3;
        double minStep = 1e-18;
        double maxStep = 600;

        double[][] tolerances = NumericalPropagator.tolerances(dPos,orbit, OrbitType.CARTESIAN);
        ODEIntegrator integrator = new DormandPrince853Integrator(minStep,maxStep,tolerances[0],tolerances[1]);
        NumericalPropagator propagator = new NumericalPropagator(integrator);

        propagator.setInitialState(initialState);
        propagator.setPositionAngleType(PositionAngle.TRUE);
        propagator.setOrbitType(OrbitType.CARTESIAN);

        // Create force models
        // NONSPHERE
        NormalizedSphericalHarmonicsProvider sphereicalHarm = GravityFieldFactory.getNormalizedProvider(degree,order);
        ForceModel nonsphere = new HolmesFeatherstoneAttractionModel(itrf,sphereicalHarm);

        // SRP
        RadiationSensitive radSens = new IsotropicRadiationSingleCoefficient(area,cR);
        SolarRadiationPressure srp = new SolarRadiationPressure(149597870700.0,1360.8/299792458 ,sun,erad,radSens);
        srp.addOccultingBody(moon,Constants.MOON_EQUATORIAL_RADIUS);

        // LUNI SOLAR
        ForceModel luni = new ThirdBodyAttraction(moon);
        ForceModel solar = new ThirdBodyAttraction(sun);

        // SOLID TIDES
        TideSystem tideSystem = sphereicalHarm.getTideSystem(); //this is the only one I could find to get this info
        ForceModel solidTide = new SolidTides(itrf,erad,earthMu,tideSystem,iers,ut1,moon,sun);

        // DRAG // DTM2000 Drag Model
        DragSensitive dragSens = new IsotropicDrag(area, dragCoeff);
        DTM2000InputParameters dtm2000params = new MarshallSolarActivityFutureEstimation(solarInd, MarshallSolarActivityFutureEstimation.StrengthLevel.AVERAGE);
        DTM2000 dtm2000 = new DTM2000(dtm2000params,sun, earth, utc);
        ForceModel drag = new DragForce(dtm2000, dragSens);

        // Add force models
        propagator.addForceModel(nonsphere);
        propagator.addForceModel(srp);
        propagator.addForceModel(luni);
        propagator.addForceModel(solar);
        propagator.addForceModel(solidTide);
        propagator.addForceModel(drag);

        return propagator;

    }

}
