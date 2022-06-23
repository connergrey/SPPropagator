import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealMatrixFormat;
import org.hipparchus.util.FastMath;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.MatricesHarvester;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitStepHandler;
import org.orekit.propagation.sampling.OrekitStepInterpolator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public class CovPropTesting {

    public static void main(String[] arg){
        // LOAD DATA AND SET CONSTANTS --------------------------------------------------------
        DataLoader loader = new DataLoader();
        loader.load(); // loads the orekit-data file to get constant parameters like tai-utc, etc.
        TimeScale utc = TimeScalesFactory.getUTC();
        String solarInd = "apr2022f10_prd.txt"; // EXAMPLE... Please update to the most revelant solar weather file
        String spVectorFileName = "/Users/connergrey/Documents/SP VECTORS/vectors_22115/scratch/SP_VEC";
        double earthMu = Constants.EGM96_EARTH_MU;
        double erad = Constants.EGM96_EARTH_EQUATORIAL_RADIUS;
        Frame eme2000 = FramesFactory.getEME2000();

        /* //INPUTS --------------------------------------------------------------------------------------
        int noradID = 50719;
        double[] firstDate = new double[]{0.51785309,392.4061074,853.4503173,6786.569731,-5.764628809,4.954299944,-0.272955645,3127.905504,-77714.88958,6460301.783,139.9940305,2259.369962,2456.54855,83.95958782,-7195.512353,-2.434589876,8.018672064,-2.996436413,55.09363219,-0.153421992,-0.058480912,0.002995817,0.062680029,-0.549127067,-1.969927331,0.000505351,-5.0957E-05,0.002485375};
        AbsoluteDate firstTCA = new AbsoluteDate("2022-04-26T06:07:33.779Z",utc);
        double[] secondDate = new double[]{0.772287095,392.3036365,853.5379877,6786.563357,-5.764641069,4.954283728,-0.273109122,4229.629853,-216558.4428,21670325.22,165.3006715,-818.3419101,2457.057498,238.6154236,-24138.21443,0.993562823,26.89143474,-3.852707506,159.7318918,-0.172874429,-0.175036557,0.003656631,0.156336716,-10.83195492,-1.967857428,0.011959412,-0.00012145,0.002492422};
*/
        //INPUTS --------------------------------------------------------------------------------------
        //Position and Velocity are in EMEJ2000
        //Covariance is in RIC
        int noradID = 50719;
        double[] firstDate = new double[]{3.702187442,352.8003634,915.248831,6782.130818,-5.480980128,5.25418048,-0.416003709,16024.98174,-5673316.01,2163463710.0,525.5483836,-130949.2328,2528.389114,6312.990729,-2408146.76,145.7972864,2680.508264,-14.5312561,5074.74461,-0.496900633,-5.646498162,0.013232742,3.943935639,-1463.532461,-0.72671436,1.628969197,-0.003533202,0.002141776};
        AbsoluteDate firstTCA = new AbsoluteDate("2022-04-29T10:33:12.267Z",utc);
        double[] secondDate = new double[]{4.541540613,371.0847151,903.515913,6782.841109,-5.425047299,5.312538698,-0.404279544,24531.86175,-11256955.91,5377397209.0,789.9911914,-292890.0894,2561.227215,12528.26547,-5985567.321,326.0464508,6662.525408,-22.86035804,10417.73311,-0.754941018,-11.59391754,0.021339026,7.656190493,-3605.527105,-0.395794437,4.013218752,-0.007092165,0.003434588};

        // [OPTIMIZE THIS!!!!] RIC PROCESS NOISE MATRIX -----------------------------------------------
        RealMatrix processNoiseRIC = MatrixUtils.createRealDiagonalMatrix(new double[]{
                4e-15,9e-15,5e-18 //0,0,0 //
        }); //in (m/s^2)^2

        //DO ALL THE STUFF --------------------------------------------------------------------------------------

        //READ THE SP VECTOR FOR THIS NORAD ID
        SPVecReader spVecReader = new SPVecReader();
        spVecReader.readSPVec(spVectorFileName,noradID);
        int degOrd = spVecReader.getDegOrd();
        double cR = spVecReader.getcR();
        double cD = spVecReader.getcD();

        //PARSE THE INPUTS AND STORE IN PROPER VARIABLE TYPES AND FRAMES -----------------------------------------

        //GET THE INITIAL POSITION AND VELOCITY AND COVARIANCE
        AbsoluteDate initialDate = firstTCA;
        Vector3D pos = new Vector3D(firstDate[1],firstDate[2],firstDate[3]).scalarMultiply(1000); //m
        Vector3D vel = new Vector3D(firstDate[4],firstDate[5],firstDate[6]).scalarMultiply(1000); //m/s
        PVCoordinates pvCoordinatesEME2000 = new PVCoordinates(pos ,vel );

        double[] first_cov_ric = new double[]{firstDate[7],firstDate[8],firstDate[9],firstDate[10],firstDate[11],firstDate[12],
                firstDate[13],firstDate[14],firstDate[15], firstDate[16],firstDate[17],firstDate[18],firstDate[19],
                firstDate[20],firstDate[21],firstDate[22],firstDate[23],firstDate[24], firstDate[25],firstDate[26],
                firstDate[27]};
        RealMatrix initCovRIC = createSymmetricMatrix(first_cov_ric);

        //DETERMINE HOW LONG TO PROPAGATE
        double firstDaysProp = firstDate[0];
        double secondDaysProp = secondDate[0];
        double propTime = (secondDaysProp - firstDaysProp) * 86400;

        //GET THE FINAL COVARIANCE GOAL
        double[] second_cov_ric = new double[]{secondDate[7],secondDate[8],secondDate[9],secondDate[10],secondDate[11],secondDate[12],
                secondDate[13],secondDate[14],secondDate[15], secondDate[16],secondDate[17],secondDate[18],secondDate[19],
                secondDate[20],secondDate[21],secondDate[22],secondDate[23],secondDate[24], secondDate[25],secondDate[26],
                secondDate[27]};
        RealMatrix goalCovRIC = createSymmetricMatrix(second_cov_ric);

        //CONVERT INITIAL COVARIANCE FROM RIC TO XYZ (EME2K) -----------------------------------------
        RealMatrix RICtoXYZ  = getRICtoXYZ(pvCoordinatesEME2000.getPosition(), pvCoordinatesEME2000.getVelocity());
        RealMatrix zeros = MatrixUtils.createRealMatrix(3,3);
        RealMatrix RICtoXYZ6D = mergeMatrices(RICtoXYZ, zeros, zeros, RICtoXYZ);
        RealMatrix initCovXYZ = RICtoXYZ6D.multiply(initCovRIC).multiplyTransposed(RICtoXYZ6D);

        double tol = 1e-5;
        while( tol < 200 ){

            //CREATE THE PROPAGATOR AND ORBIT -----------------------------------------------------------------

            //CREATE ORBIT
            CartesianOrbit orbit = new CartesianOrbit(pvCoordinatesEME2000, eme2000, initialDate, earthMu);
            //CREATE THE PROPAGATOR
            NumericalPropagator propagator = PropCreator.createProp(orbit, earthMu, erad, utc, degOrd, cR, cD, solarInd, tol);
            //SET MATRICES HARVESTER
            MatricesHarvester harvester = propagator.setupMatricesComputation("stm", null, null);
            //SET STEP HANDLER
            StepHandler stepHandler = new StepHandler(initCovXYZ, processNoiseRIC, harvester);
            propagator.setStepHandler(stepHandler);

            //PROPAGATE -----------------------------------------------------------------
            //PROPAGATE STATE
            SpacecraftState state = propagator.propagate(initialDate.shiftedBy(propTime));

            //GET PROPAGATED COVARIANCE MATRIX
            RealMatrix finalCovXYZ = stepHandler.getCov();

            //TRANSFORM THE FINAL COVARIANCE TO THE RIC FRAME

            RealMatrix finalRICtoXYZ = getRICtoXYZ(state.getPVCoordinates().getPosition(), state.getPVCoordinates().getVelocity());
            RealMatrix finalXYZtoRIC = finalRICtoXYZ.transpose();
            RealMatrix XYZtoRIC6D = mergeMatrices(finalXYZtoRIC, zeros, zeros, finalXYZtoRIC);
            RealMatrix finalCovRIC = XYZtoRIC6D.multiply(finalCovXYZ).multiplyTransposed(XYZtoRIC6D);

            // OUTPUTS
            RealMatrixFormat matrixFormat = new RealMatrixFormat("", "", "", "\n", "", ", ");

        /*System.out.println("Cost Function:");
        System.out.println(finalCovRIC.getSubMatrix(0, 2, 0, 2).subtract(goalCovRIC.getSubMatrix(0, 2, 0, 2)).getFrobeniusNorm());
        System.out.println();

        System.out.println("PROPAGATED MATRIX: [m^2]");
        System.out.println(matrixFormat.format(finalCovRIC.getSubMatrix(0,2,0,2)));
        System.out.println( "POSITION STD: " + FastMath.sqrt( finalCovRIC.getSubMatrix(0,2,0,2).getTrace()) + " meters");
        System.out.println();
        System.out.println();

        System.out.println("GOAL MATRIX: [m^2]");
        System.out.println(matrixFormat.format(goalCovRIC.getSubMatrix(0,2,0,2)));
        System.out.println( "POSITION STD: " + FastMath.sqrt( goalCovRIC.getSubMatrix(0,2,0,2).getTrace()) + " meters");*/

            System.out.println("TOLERANCE: " + tol + " -> " +
                    "POSITION STD: " + FastMath.sqrt(finalCovRIC.getSubMatrix(0, 2, 0, 2).getTrace()) + " meters");
            System.out.println("POSITION COVARIANCE MATRIX:");
            System.out.println(matrixFormat.format(finalCovRIC.getSubMatrix(0,2,0,2)));
            System.out.println();

            tol = tol*10;
        }
    /*    // Adding the process noise at the end
        // compute the stm from previous time to current time
        RealMatrix stm = harvester.getStateTransitionMatrix(state);
        // propagate the covariance
        RealMatrix propCovXYZ = stm.multiply(initCovXYZ).multiplyTransposed(stm);
        // transform to RIC
        RealMatrix propCovRIC = XYZtoRIC6D.multiply(propCovXYZ).multiplyTransposed(XYZtoRIC6D);
        // add the process noise
        RealMatrix PNtrans = StepHandler.processNoiseTransition(propTime);
        // add
        RealMatrix covWNoiseRIC = propCovRIC.add( PNtrans.multiply(processNoiseRIC).multiplyTransposed(PNtrans) );

        System.out.println("JUST PROP: [m^2]");
        System.out.println(matrixFormat.format(propCovRIC.getSubMatrix(0,2,0,2)));
        System.out.println( "POSITION STD: " + FastMath.sqrt( propCovRIC.getSubMatrix(0,2,0,2).getTrace()) + " meters");

        System.out.println("PROP W NOISE: [m^2]");
        System.out.println(matrixFormat.format(covWNoiseRIC.getSubMatrix(0,2,0,2)));
        System.out.println( "POSITION STD: " + FastMath.sqrt( covWNoiseRIC.getSubMatrix(0,2,0,2).getTrace()) + " meters");

        System.out.println("JUST NOISE: [m^2]");
        System.out.println(matrixFormat.format(PNtrans.multiply(processNoiseRIC).multiplyTransposed(PNtrans).getSubMatrix(0,2,0,2)));
        System.out.println( "POSITION STD: " + FastMath.sqrt( PNtrans.multiply(processNoiseRIC).multiplyTransposed(PNtrans).getSubMatrix(0,2,0,2).getTrace()) + " meters");
    */
    }

    public static class StepHandler implements OrekitStepHandler {

        private RealMatrix cov;
        private RealMatrix stm0;
        private RealMatrix stm2;
        private RealMatrix stm;
        private final RealMatrix processNoiseRIC;
        private AbsoluteDate lastDate;
        private MatricesHarvester harvester;

        public StepHandler(RealMatrix initialCovariance, RealMatrix processNoiseRIC, MatricesHarvester harvester){
            this.cov = initialCovariance;
            this.processNoiseRIC = processNoiseRIC;
            this.harvester = harvester;

        }

        @Override
        public void init(final SpacecraftState s0, final AbsoluteDate t) {
            stm0 = MatrixUtils.createRealIdentityMatrix(6);
            lastDate = t;
        }

        @Override
        public void handleStep(OrekitStepInterpolator interpolator) {
            // get the state at current time
            SpacecraftState state = interpolator.getCurrentState();
            // compute the stm from previous time to current time
            stm2 = harvester.getStateTransitionMatrix(state);
            // calculate the stm via stm chain rule
            stm = stm2.multiply( MatrixUtils.inverse(stm0) );
            // propagate the covariance
            cov = stm.multiply(cov).multiplyTransposed(stm);

            // add process noise
            // find the RIC to XYZ transformation
            Vector3D pos = state.getPVCoordinates().getPosition();
            Vector3D vel = state.getPVCoordinates().getVelocity();
            RealMatrix RICtoXYZ = getRICtoXYZ(pos,vel);

            //transform the process noise to XYZ
            RealMatrix processNoiseXYZ = RICtoXYZ.multiply(processNoiseRIC).multiplyTransposed(RICtoXYZ);
            // determine how long we propagated
            double dt = state.getDate().durationFrom( lastDate );
            // compute gamma
            RealMatrix gamma = processNoiseTransition(dt);
            // compute process noise covariance
            RealMatrix processNoiseCovarianceXYZ = gamma.multiply(processNoiseXYZ).multiplyTransposed(gamma);
            // add the process noise
            cov = cov.add( processNoiseCovarianceXYZ );

            // update stm0 and last date
            stm0 = stm2;
            lastDate = state.getDate();

        }

        @Override
        public void finish(SpacecraftState finalState) {


        }

        public RealMatrix getCov(){
            return cov;
        }

        public static RealMatrix processNoiseTransition(double dt){

            RealMatrix gamma = MatrixUtils.createRealMatrix(6,3);

            double d1 = dt*dt/2;
            double d2 = dt;

            // set the upper diagonal
            gamma.setEntry(0,0,d1);
            gamma.setEntry(1,1,d1);
            gamma.setEntry(2,2,d1);

            // set the lower diagonal
            gamma.setEntry(3,0,d2);
            gamma.setEntry(4,1,d2);
            gamma.setEntry(5,2,d2);

            return gamma;
        }


    }

    public static RealMatrix getRICtoXYZ(Vector3D pos, Vector3D vel){

        // define the RIC coordinate directions
        // define screen's RIC frame (same as RTN)
        Vector3D rHat = pos.normalize();
        Vector3D hVec = pos.crossProduct(vel);
        Vector3D cHat = hVec.normalize();
        Vector3D iHat = cHat.crossProduct(rHat);

        // create RIC to XYZ transformation matrix
        RealMatrix XYZtoRIC = MatrixUtils.createRealMatrix( new double[][] { {rHat.getX(), rHat.getY(), rHat.getZ()},
                {iHat.getX(), iHat.getY(), iHat.getZ()} , {cHat.getX(), cHat.getY(), cHat.getZ()} } );

        RealMatrix RICtoXYZ = XYZtoRIC.transpose();

        return RICtoXYZ;

    }

    public static RealMatrix createSymmetricMatrix(double[] mat){

        // input mat in order: RR	IR	II	CR	CI	CC
        //            CRDOT_R
        //            CRDOT_T
        //            CRDOT_N
        //            CRDOT_RDOT
        //            CTDOT_R
        //            CTDOT_T
        //            CTDOT_N
        //            CTDOT_RDOT
        //            CTDOT_TDOT
        //            CNDOT_R
        //            CNDOT_T
        //            CNDOT_N
        //            CNDOT_RDOT
        //            CNDOT_TDOT
        //            CNDOT_NDOT

        double RR = mat[0];
        double TR = mat[1];
        double TT = mat[2];
        double NR = mat[3];
        double NT = mat[4];
        double NN = mat[5];

        double RDOT_R = mat[6];
        double RDOT_T = mat[7];
        double RDOT_N = mat[8];
        double RDOT_RDOT = mat[9];
        double TDOT_R = mat[10];
        double TDOT_T = mat[11];
        double TDOT_N = mat[12];
        double TDOT_RDOT = mat[13];
        double TDOT_TDOT = mat[14];
        double NDOT_R = mat[15];
        double NDOT_T = mat[16];
        double NDOT_N = mat[17];
        double NDOT_RDOT = mat[18];
        double NDOT_TDOT = mat[19];
        double NDOT_NDOT = mat[20];

        //output matrix: [ [RR, TR, NR]               [RDOT_R, TDOT_R, NDOT_R]
        //                 [TR, TT, NT]               [RDOT_T, TDOT_T, NDOT_T]
        //                 [NR, NT, NN]               [RDOT_N, TDOT_N, NDOT_N]

        //                 [RDOT_R, RDOT_T, RDOT_N]   [RDOT_RDOT, TDOT_RDOT, NDOT_RDOT]
        //                 [TDOT_R, TDOT_T, TDOT_N]   [TDOT_RDOT, TDOT_TDOT, NDOT_TDOT]
        //                 [NDOT_R, NDOT_T, NDOT_N]   [NDOT_RDOT, NDOT_TDOT, NDOT_NDOT] ]
        return MatrixUtils.createRealMatrix(new double[][]
                {
                        {RR, TR, NR, RDOT_R, TDOT_R, NDOT_R},
                        {TR, TT, NT, RDOT_T, TDOT_T, NDOT_T},
                        {NR, NT, NN, RDOT_N, TDOT_N, NDOT_N},
                        {RDOT_R, RDOT_T, RDOT_N, RDOT_RDOT, TDOT_RDOT, NDOT_RDOT},
                        {TDOT_R, TDOT_T, TDOT_N, TDOT_RDOT, TDOT_TDOT, NDOT_TDOT},
                        {NDOT_R, NDOT_T, NDOT_N, NDOT_RDOT, NDOT_TDOT, NDOT_NDOT}
                });

    }

    public static RealMatrix mergeMatrices(RealMatrix TLxyz,RealMatrix TRxyz,RealMatrix BLxyz,RealMatrix BRxyz){

        RealMatrix initCovXYZ = MatrixUtils.createRealMatrix(6,6);
        for (int k = 1; k < 5; k++) {
            RealMatrix mat = null;
            int i0 = 0;
            int j0 = 0;
            if(k==1){
                mat = TLxyz;
            }
            else if(k==2){
                mat = TRxyz;
                i0 = 0;
                j0 = 3;
            }
            else if(k==3){
                mat = BLxyz;
                i0 = 3;
                j0 = 0;}
            else if(k==4){
                mat = BRxyz;
                i0 = 3;
                j0 = 3;}

            for (int i = 0; i < TLxyz.getData().length; i++) {
                for (int j = 0; j < TLxyz.getData().length; j++) {

                    initCovXYZ.setEntry(i0+i,j0+j, mat.getEntry(i,j));
                }
            }
        }
        return initCovXYZ;
    }


}
