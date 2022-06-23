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

import java.io.FileNotFoundException;

public class stmTesting {

    public static void main(String[] arg) throws FileNotFoundException {

        // make sure you update orekit-data

        DataLoader loader = new DataLoader();
        loader.load(); // loads the orekit-data file to get constant parameters like tai-utc, etc.

        stmWork stmwork = new stmWork();
        stmwork.work();


        /*String solarInd = "mar2022f10_prd.txt"; // EXAMPLE... Please update to the most revelant solar weather file
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
        NumericalPropagator propagator = PropCreator.createProp(orbit, earthMu, erad, utc, degOrd, cR, cD, solarInd);

        // specify the initial covariance
        RealMatrix initialCovariance = MatrixUtils.createRealIdentityMatrix(6);
        initialCovariance.setEntry(3,3,0.001);
        initialCovariance.setEntry(4,4,0.001); // set the velocity uncertainty lower
        initialCovariance.setEntry(5,5,0.001);

        // set up matricies harvester
        MatricesHarvester harvester = propagator.setupMatricesComputation("stm",null, null);

        // propagate
        double propTime = 20;
        SpacecraftState state = propagator.propagate(initialDate.shiftedBy(propTime));


        RealMatrix stm = harvester.getStateTransitionMatrix( state );
        RealMatrix finalCovariance = stm.multiply(initialCovariance).multiplyTransposed(stm);
        System.out.println(finalCovariance);*/

/*        RealMatrixFormat matrixFormat = new RealMatrixFormat("", "", "", "\n", "", ", ");

        System.out.println(matrixFormat.format(finalCovariance));
        System.out.println(matrixFormat.format(stm));

        */
/*

        ///////////// test if the STM works ////////////////

        //Vector3D finalPos = state.getPVCoordinates().getPosition();
        //Vector3D finalVel = state.getPVCoordinates().getVelocity();

        //1. deviate the initial orbit and see how it impacts the final

        //deviation from the initial orbit
        Vector3D delx = new Vector3D(25,25,-25);
        Vector3D delv = new Vector3D(0,0,0);

        //modify the initial condition
        //Vector3D initPrime = propagator.getInitialState().getPVCoordinates().getPosition().add(delx);
        PVCoordinates pvCoordinatesEME2000Prime = new PVCoordinates(initial.getPosition().add(delx), initial.getVelocity().add(delv));
        CartesianOrbit orbitPrime = new CartesianOrbit(pvCoordinatesEME2000Prime, eme2000, initialDate, earthMu);
        NumericalPropagator propagatorPrime = PropCreator.createProp(orbitPrime, earthMu, erad, utc, degOrd, cR, cD, solarInd);

        //propagate the same amount of time
        SpacecraftState statePrime = propagatorPrime.propagate(initialDate.shiftedBy(propTime));
        Vector3D finalPosPrime =  statePrime.getPVCoordinates().getPosition();
        Vector3D finalDelPos = finalPosPrime.subtract( finalPos );
        Vector3D finalVelPrime =  statePrime.getPVCoordinates().getVelocity();
        Vector3D finalDelVel = finalVelPrime.subtract( finalVel );

        //2. use the stm to compute the final deviation

        //create deviation column vector
        RealMatrix stateDeviation = MatrixUtils.createColumnRealMatrix(new double[]{delx.getX(),delx.getY(),delx.getZ(),delv.getX(),delv.getY(),delv.getZ()});
        RealMatrix stmFinalDev = stm.multiply(stateDeviation);

        //3. Compare
        //System.out.println(finalDelPos + "," + finalDelVel);
        //System.out.println( stmFinalDev );
*/

    }

    public static class stmWork{

        public static MatricesHarvester harvester;
        public static NumericalPropagator propagator;

        public void work(){

            String solarInd = "mar2022f10_prd.txt"; // EXAMPLE... Please update to the most revelant solar weather file
            String spVectorFileName = "/Users/connergrey/Documents/SP VECTORS/vectors_22115/scratch/SP_VEC/50/50719"; // EXAMPLE... SP VECTOR FILE NAME

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

            // get the state and time at TCA
            AbsoluteDate initialDate = new AbsoluteDate("2022-04-29T10:33:12.267Z",utc);
            Vector3D pos = new Vector3D(352.8003634,	915.248831,	6782.130818).scalarMultiply(1000);
            Vector3D vel = new Vector3D(-5.480980128,	5.25418048,	-0.416003709).scalarMultiply(1000);
            PVCoordinates pvCoordinatesEME2000 = new PVCoordinates(pos ,vel );

            // specify initial covariance in RIC
            double[] inp_mat = new double[]{16024.98174,-5673316.01,2163463710.0,525.5483836,-130949.2328,2528.389114,
                    6312.990729,-2408146.76,145.7972864,2680.508264,-14.5312561,5074.74461,-0.496900633,-5.646498162,
                    0.013232742,3.943935639,-1463.532461,-0.72671436,1.628969197,-0.003533202,0.002141776};
            RealMatrix initCovRIC = createSymmetricMatrix(inp_mat);

            // process noise in RIC
            RealMatrix processNoiseRIC = MatrixUtils.createRealDiagonalMatrix(new double[]{
                    4e-15,9e-15,5e-18
            }); //in (m/s^2)^2

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

            // transform covariance to XYZ
            // need to create 3 sub matricies and transform them, then rejoin somehow
            RealMatrix TL = initCovRIC.getSubMatrix(0,2,0,2);
            RealMatrix TR = initCovRIC.getSubMatrix(0,2,3,5);
            RealMatrix BR = initCovRIC.getSubMatrix(3,5,3,5);

            RealMatrix TLxyz = RICtoXYZ.multiply(TL).multiplyTransposed(RICtoXYZ);// k=1
            RealMatrix TRxyz = RICtoXYZ.multiply(TR).multiplyTransposed(RICtoXYZ);// k=2
            RealMatrix BLxyz = TRxyz.transpose();// k=3
            RealMatrix BRxyz = RICtoXYZ.multiply(BR).multiplyTransposed(RICtoXYZ);// k=4

            // merge them into the 6x6 matrix
            RealMatrix initCovXYZ = mergeMatrices(TLxyz, TRxyz, BLxyz, BRxyz);

            RealMatrixFormat matrixFormat = new RealMatrixFormat("", "", "", "\n", "", ", ");
            //System.out.println(matrixFormat.format(initCovXYZ.getSubMatrix(3,5,3,5)));
            //System.out.println(matrixFormat.format(initCovXYZ.getSubMatrix(0,2,0,2)));
            //System.out.println( FastMath.sqrt( initCovXYZ.getTrace()));

            // Create orbit from initial condition
            CartesianOrbit orbit = new CartesianOrbit(pvCoordinatesEME2000, eme2000, initialDate, earthMu);

            // Create the integrator and propagator with desired models
            propagator = PropCreator.createProp(orbit, earthMu, erad, utc, degOrd, cR, cD, solarInd);


            // set up matricies harvester
            harvester = propagator.setupMatricesComputation("stm",null, null);

            // set up step handler
            propagator.setStepHandler( new StepHandler(initCovXYZ, processNoiseRIC) );

            // propagate
            double propDays = 4.541540613 - 3.702187442;
            double propTime = propDays * 86400;
            //double propTime = 1000;
            SpacecraftState state = propagator.propagate(initialDate.shiftedBy(propTime));

            // propagate covariance
            //RealMatrix stm = harvester.getStateTransitionMatrix(state);
            //RealMatrix finalCov = stm.multiply(initCovXYZ).multiplyTransposed(stm);

            // grab the propagated covariance
            RealMatrix finalCovXYZ = StepHandler.finalCovXYZ;

            // transform to RIC
                // define the RIC coordinate directions
                Vector3D pos1 = state.getPVCoordinates().getPosition();
                Vector3D vel1 = state.getPVCoordinates().getVelocity();
                // define screen's RIC frame (same as RTN)
                Vector3D rHat1 = pos1.normalize();
                Vector3D hVec1 = pos1.crossProduct(vel1);
                Vector3D cHat1 = hVec1.normalize();
                Vector3D iHat1 = cHat1.crossProduct(rHat1);

                // create RIC to XYZ transformation matrix
                RealMatrix XYZtoRICFinal = MatrixUtils.createRealMatrix( new double[][] { {rHat1.getX(), rHat1.getY(), rHat1.getZ()},
                        {iHat1.getX(), iHat1.getY(), iHat1.getZ()} , {cHat1.getX(), cHat1.getY(), cHat1.getZ()} } );

                // need to create 3 sub matricies and transform them, then rejoin somehow
                RealMatrix TLf = finalCovXYZ.getSubMatrix(0,2,0,2);
                RealMatrix TRf = finalCovXYZ.getSubMatrix(0,2,3,5);
                RealMatrix BRf = finalCovXYZ.getSubMatrix(3,5,3,5);

                RealMatrix TLfric = XYZtoRICFinal.multiply( TLf ).multiplyTransposed(XYZtoRICFinal);// k=1
                RealMatrix TRfric = XYZtoRICFinal.multiply( TRf ).multiplyTransposed(XYZtoRICFinal);// k=2
                RealMatrix BLfric = TRxyz.transpose();// k=3
                RealMatrix BRfric = XYZtoRICFinal.multiply( BRf ).multiplyTransposed(XYZtoRICFinal);// k=4

                // merge them into the 6x6 matrix
                RealMatrix finalCovRIC = mergeMatrices(TLfric, TRfric, BLfric, BRfric);


            System.out.println("MINE:");
            System.out.println(matrixFormat.format(finalCovRIC.getSubMatrix(0,2,0,2)));
            System.out.println( FastMath.sqrt( finalCovRIC.getTrace()));


            // goal covariance I am shooting for
            // specify initial covariance in RIC
            double[] goal_mat = new double[]{24531.86175,-11256955.91,5377397209.0,789.9911914,-292890.0894,
                    2561.227215,12528.26547,-5985567.321,326.0464508,6662.525408,-22.86035804,10417.73311,
                    -0.754941018,-11.59391754,0.021339026,7.656190493,-3605.527105,-0.395794437,4.013218752,
                    -0.007092165,0.003434588};
            RealMatrix goalFinalCovRIC = createSymmetricMatrix(goal_mat);
            System.out.println("ACTUAL:");
            System.out.println(matrixFormat.format(goalFinalCovRIC.getSubMatrix(0,2,0,2)));
            System.out.println( FastMath.sqrt( goalFinalCovRIC.getTrace()));



        }

        public RealMatrix mergeMatrices(RealMatrix TLxyz,RealMatrix TRxyz,RealMatrix BLxyz,RealMatrix BRxyz){

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

        public RealMatrix createSymmetricMatrix(double[] mat){

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


    }



    public static class StepHandler implements OrekitStepHandler {

        private RealMatrix cov;
        private RealMatrix stm0;
        private RealMatrix stm2;
        private RealMatrix stm;
        private RealMatrix processNoiseRIC;
        private AbsoluteDate lastDate;
        public static RealMatrix finalCovXYZ;

        public StepHandler(RealMatrix initialCovariance,RealMatrix processNoiseRIC){
            this.cov = initialCovariance;
            this.processNoiseRIC = processNoiseRIC;

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
            stm2 = stmWork.harvester.getStateTransitionMatrix(state);

            // calculate the stm via stm chain rule
            stm = stm2.multiply( MatrixUtils.inverse(stm0) );

            // propagate the covariance
            this.cov = stm.multiply(this.cov).multiplyTransposed(stm);

            // add process noise

                // find the RIC to XYZ transformation

                    Vector3D pos = state.getPVCoordinates().getPosition();
                    Vector3D vel = state.getPVCoordinates().getVelocity();
                    // define screen's RIC frame (same as RTN)
                    Vector3D rHat = pos.normalize();
                    Vector3D hVec = pos.crossProduct(vel);
                    Vector3D cHat = hVec.normalize();
                    Vector3D iHat = cHat.crossProduct(rHat);

                    // create RIC to XYZ transformation matrix
                    RealMatrix XYZtoRIC = MatrixUtils.createRealMatrix( new double[][] { {rHat.getX(), rHat.getY(), rHat.getZ()},
                            {iHat.getX(), iHat.getY(), iHat.getZ()} , {cHat.getX(), cHat.getY(), cHat.getZ()} } );

                RealMatrix RICtoXYZ = XYZtoRIC.transpose();

                //transform the process noise to XYZ
                RealMatrix processNoiseXYZ = RICtoXYZ.multiply(processNoiseRIC).multiplyTransposed(RICtoXYZ);

                // determine how long we propagated
                double dt = state.getDate().durationFrom( lastDate );

                // compute gamma
                RealMatrix gamma = processNoiseTransition(dt);

                // compute process noise covariance
                RealMatrix processNoiseCovarianceXYZ = gamma.multiply(processNoiseXYZ).multiplyTransposed(gamma);

                // add the process noise
            this.cov = this.cov.add( processNoiseCovarianceXYZ );


            // update stm0 and last date
            stm0 = stm2;
            lastDate = state.getDate();

        }

        @Override
        public void finish(SpacecraftState finalState) {

/*            RealMatrixFormat matrixFormat = new RealMatrixFormat("", "", "", "\n", "", ", ");
            System.out.println(matrixFormat.format(this.cov.getSubMatrix(0,2,0,2)));
            System.out.println( FastMath.sqrt( this.cov.getTrace()));*/
            finalCovXYZ = this.cov;

        }

        public RealMatrix processNoiseTransition(double dt){

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


}