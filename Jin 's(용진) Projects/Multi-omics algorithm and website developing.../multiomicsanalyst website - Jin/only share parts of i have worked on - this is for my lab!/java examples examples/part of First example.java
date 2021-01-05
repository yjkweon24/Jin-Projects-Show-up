import general.utils.DataUtils;
import general.utils.RDataUtils;
import org.rosuda.REngine.Rserve.RConnection;


public class ProcessBean implements Serializable {

    private final SessionBean1 sb = (SessionBean1) DataUtils.findBean("sessionBean1");  
    
    private int minCutoff = 5;


    public int getMinCutoff() {
        return minCutoff;
    }

    public void setMinCutoff(int minCutoff) {
        this.minCutoff = minCutoff;
    }

    private boolean normPerformed = false;

    public boolean isNormPerformed() {
        return normPerformed;
    }

    private String lowFilterOpt = "prevalence";
    private String varFilterOpt = "iqr";

    
    
    //low count filter 
    public String getLowFilterOpt() {
        return lowFilterOpt;
    }

    public void setLowFilterOpt(String lowFilterOpt) {
        this.lowFilterOpt = lowFilterOpt;
    }

    
    
    
    //low variance filter 
    public String getVarFilterOpt() {
        return varFilterOpt;
    }

    public void setVarFilterOpt(String varFilterOpt) {
        this.varFilterOpt = varFilterOpt;
    }

    //low count filter min count
    private int minCount = 4;

    public int getMinCount() {
        return minCount;
    }

    public void setMinCount(int minCount) {
        this.minCount = minCount;
    }

    
    private int smplPerct = 20;
    private int featPerct = 10;

    //low count filter prevalence in sample
    public int getSmplPerct() {
        return smplPerct;
    }

    public void setSmplPerct(int smplPerct) {
        this.smplPerct = smplPerct;
    }
    
    
    //normalization generalized log transformation
    private float lambda1 = 0.5f;

    public float getLambda1() {
        return lambda1;
    }
 
    public void setLambda1(float lambda1) {
        this.lambda1 = lambda1;
    }
    
    //normalization twoparameter boxcox transformation
    private float lambda1b = 0.5f;

    public float getLambda1b() {
        return lambda1b;
    }
 
    public void setLambda1b(float lambda1b) {
        this.lambda1b = lambda1b;
    }
    
    private float lambda2b = 0.5f;

    public float getLambda2b() {
        return lambda2b;
    }
 
    public void setLambda2b(float lambda2b) {
        this.lambda2b = lambda2b;
    }
    
    
    //# of variables for pick one in o2pls page for the update
    private int joint = 2;

    public int getJoint() {
        return joint;
    }
 
    public void setJoint(int joint) {
        this.joint = joint;
    }
        
    private int orth1 = 2;

    public int getOrth1() {
        return orth1;
    }
 
    public void setOrth1(int orth1) {
        this.orth1 = orth1;
    }  
    
    private int orth2 = 2;

    public int getOrth2() {
        return orth2;
    }
 
    public void setOrth2(int orth2) {
        this.orth2 = orth2;
    }    
    
    private int jointc = 6;

    public int getJointc() {
        return jointc;
    }
 
    public void setJointc(int jointc) {
        this.jointc = jointc;
    }
        
    private int orth1c = 8;

    public int getOrth1c() {
        return orth1c;
    }
 
    public void setOrth1c(int orth1c) {
        this.orth1c = orth1c;
    }  
    
    private int orth2c = 8;

    public int getOrth2c() {
        return orth2c;
    }
 
    public void setOrth2c(int orth2c) {
        this.orth2c = orth2c;
    }        


    //loop # choose for fixing joint and orth
    private int loops = 20;

    public int getLoops() {
        return loops;
    }
 
    public void setLoops(int loops) {
        this.loops = loops;
    }        
    
    
    
    //# of variables for pick one in o2pls page for the submit
    public int jointb = 2;

    public int getJointb() {
        return jointb;
    }
 
    public void setJointb(int jointb) {
        this.jointb = jointb;
    }
        
    public int orth1b = 2;

    public int getOrth1b() {
        return orth1b;
    }
 
    public void setOrth1b(int orth1b) {
        this.orth1b = orth1b;
    }  
    
    public int orth2b = 2;

    public int getOrth2b() {
        return orth2b;
    }
 
    public void setOrth2b(int orth2b) {
        this.orth2b = orth2b;
    }    
       
    //PA analysis 
    private int perm = 999;

    public int getPerm() {
        return perm;
    }
 
    public void setPerm(int perm) {
        this.perm = perm;
    }
    
    
    
    
    private float size = 0.5f;

    public float getSize() {
        return size;
    }
 
    public void setSize(float size) {
        this.size = size;
    }
    
    
    
    private int type = 1;
    
    
    public int getType() {
        return type;
    }

    public void setType(int type) {
        this.type = type;
    }    
    
    
    private int color = 1;
    
    
    public int getColor() {
        return color;
    }

    public void setColor(int color) {
        this.color = color;
    }  
    
    
    private int arrow = 1;
    
    
    public int getArrow() {
        return arrow;
    }

    public void setArrow(int arrow) {
        this.arrow = arrow;
    }  
    
    
    
    //PA Update
    private int choice = 1;
    
    
    public int getChoice() {
        return choice;
    }

    public void setChoice(int choice) {
        this.choice = choice;
    }  
    
    
    private int adjustment = 1;
    
    
    public int getAdjustment() {
        return adjustment;
    }

    public void setAdjustment(int adjustment) {
        this.adjustment = adjustment;
    }    
    
    
    
    private int perm2 = 999;

    public int getPerm2() {
        return perm2;
    }
 
    public void setPerm2(int perm2) {
        this.perm2 = perm2;
    }
    
       
       
    
    
    //low variance filter % to remove
    public int getFeatPerct() {
        return featPerct;
    }

    public void setFeatPerct(int featPerct) {
        this.featPerct = featPerct;
    }

    private double cov = 3.0;

    public double getCov() {
        return cov;
    }

    public void setCov(double cov) {
        this.cov = cov;
    }

    
    //used in ... disabled ="#{!procBean.annotated}"
    private boolean annotated = false;

    public boolean isAnnotated() {
        return annotated;
    }

    public void setAnnotated(boolean annotated) {
        this.annotated = annotated;
    }

    public String prepareDataFilter() {
        RConnection RC = sb.getRConnection();
        String format = sb.getDataFormat();
        String datatype = sb.getModuleType();
        
        int res = RDataUtils.missingCheckData(RC, datatype);
        
        String id = sb.getModuleType();
        if (id.equals("omics")) {


            return "Data Filter1";
        }
        FacesContext.getCurrentInstance().addMessage(null,
                new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to proceed!"));
        return null;
    }
    
    
    //this is for the omcis1 data missing value imputations
    public void performMissing1() {
        annotated = false;
        RConnection RC = sb.getRConnection();
        String missImp = sb.getMissImp();
        
        
        int res = RDataUtils.applyMissing(RC, missImp);
        
        if (res == 1) {
            annotated = true;
            String idOptMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Missing value imputation - OK", idOptMsg));
            sb.resetPageInit("Missing value1");
        } else {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform Missing value imputations!"));
        }
    }  
    
    //this is for the omcis2 data missing value imputations
    public void performMissing2() {
        annotated = false;
        RConnection RC = sb.getRConnection();
        String missImp = sb.getMissImp();
        
        
        int res = RDataUtils.applyMissing2(RC, missImp);
        
        if (res == 1) {
            annotated = true;
            String idOptMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Missing value imputation - OK", idOptMsg));
            sb.resetPageInit("Missing value2");
        } else {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform Missing value imputations!"));
        }
    } 
    
    //this is for the both omcis1 and 2 data missing value imputations
    public void performMissing3() {
        annotated = false;
        RConnection RC = sb.getRConnection();
        String missImp = sb.getMissImp();
        String missImp2 = sb.getMissImp2();
        
        int res = RDataUtils.applyMissing3(RC, missImp, missImp2);
        
        if (res == 1) {
            annotated = true;
            String idOptMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Missing value imputation - OK", idOptMsg));
            sb.resetPageInit("Missing value12");
        } else {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform Missing value imputations!"));
        }
    }         
    

    //For the first omics data filter
    public void performFilteration() {
        annotated = false;
        RConnection RC = sb.getRConnection();
        int res = 0;
        
        //low count
        res = RDataUtils.applyAbundanceFilter(RC, lowFilterOpt, minCount, smplPerct / 100.0);
        if (res == 0) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform low abundance filtering!"));
            return;
        }
        
        //low variance 
        res = RDataUtils.applyVarianceFilter(RC, varFilterOpt, featPerct / 100.0);
        if (res == 0) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform low variance filtering!"));
            return;
        }

        String idOptMsg = RDataUtils.getCurrentMsg(RC);
        FacesContext.getCurrentInstance().addMessage(null,
                new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Filtering- OK", idOptMsg));
        annotated = true;

        sb.resetPageInit("Data Filter1");
    }
    
    

    //For the second omics data filter
    public void performFilteration2() {
        annotated = false;
        RConnection RC = sb.getRConnection();
        int res = 0;
        
        //low count
        res = RDataUtils.applyAbundanceFilter2(RC, lowFilterOpt, minCount, smplPerct / 100.0);
        if (res == 0) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform low abundance filtering!"));
            return;
        }
        
        //low variance 
        res = RDataUtils.applyVarianceFilter2(RC, varFilterOpt, featPerct / 100.0);
        if (res == 0) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform low variance filtering!"));
            return;
        }

        String idOptMsg = RDataUtils.getCurrentMsg(RC);
        FacesContext.getCurrentInstance().addMessage(null,
                new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Filtering- OK", idOptMsg));
        annotated = true;

        sb.resetPageInit("Data Filter2");
    }    
    
    
    //this is for the first omics normalization.
    public void performNormalization() {
        RConnection RC = sb.getRConnection();
        //String rareOpt = sb.getRareOpt();
        String scaleOpt = sb.getScaleOpt();
        String transformOpt = sb.getTransformOpt();
        
        if (!scaleOpt.equals("none") && !transformOpt.equals("none")) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Cannot perform data scaling or data transformation together! Please choose one of them."));
        }
        
        int res = RDataUtils.performNormalization(RC, scaleOpt, transformOpt, lambda1, lambda1b, lambda2b);
        if (res == 1) {
            normPerformed = true;
            String idOptMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Normalization - OK", idOptMsg));
            sb.resetPageInit("Normalization1");
        } else {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform normalization!"));
        }
    }
    
    //this is for the second omics normalization.
    public void performNormalization2() {
        RConnection RC = sb.getRConnection();
        //String rareOpt = sb.getRareOpt();
        String scaleOpt = sb.getScaleOpt();
        String transformOpt = sb.getTransformOpt();
        
        if (!scaleOpt.equals("none") && !transformOpt.equals("none")) {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Cannot perform data scaling or data transformation together! Please choose one of them."));
        }
        
        int res = RDataUtils.performNormalization2(RC, scaleOpt, transformOpt, lambda1, lambda1b, lambda2b);
        if (res == 1) {
            normPerformed = true;
            String idOptMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_INFO, "Data Normalization - OK", idOptMsg));
            
            int[] res2 = RDataUtils.summarydata(RC);
            
            String msg = "<table>";
            msg = msg + "<tr> <td align=\"left\" width=\"250\"><td align=\"left\"><b>Sample</b></td> <td align=\"center\" width=\"70\"><b>Feature</b></td> </tr>";
            msg = msg + "<tr> <td align=\"left\"> <b>=================</b></td> </tr>";
            msg = msg + "<tr><td align=\"left\"> <b>Omics1- original:</b></td><td align=\"center\"> " + res2[1] + "</td><td align=\"center\"width=\"60\">" + res2[3] + "</td></tr>";
            msg = msg + "<tr><td align=\"left\"> <b>Omics2- original:</b></td><td align=\"center\"> " + res2[2] + "</td><td align=\"center\"width=\"60\">" + res2[4] + "</td></tr>";
            msg = msg + "<tr> <td align=\"left\"> <b>=================</b></td> </tr>";
            msg = msg + "<tr><td align=\"left\"> <b>Omics1- Filtered:</b></td><td align=\"center\"> " + res2[5] + "</td><td align=\"center\"width=\"60\">" + res2[7] + "</td></tr>";
            msg = msg + "<tr><td align=\"left\"> <b>Omics2- Filtered:</b></td><td align=\"center\"> " + res2[6] + "</td><td align=\"center\"width=\"60\">" + res2[8] + "</td></tr>";
            msg = msg + "<tr> <td align=\"left\"> <b>=================</b></td> </tr>";
            msg = msg + "<tr><td align=\"left\"> <b>Omics1- Normalized:</b></td><td align=\"center\"> " + res2[9] + "</td><td align=\"center\"width=\"60\">" + res2[11] + "</td></tr>";
            msg = msg + "<tr><td align=\"left\"> <b>Omics2- Normalized:</b></td><td align=\"center\"> " + res2[10] + "</td><td align=\"center\"width=\"60\">" + res2[12] + "</td></tr>";
            msg = msg + "</table>";

            sb.setProcMsg2(msg); //set table procmsg2 

            sb.resetPageInit("Normalization2"); // return to the page
        } else {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform normalization!"));
        }
    }
    
    //Each has two: for two omics.... This comes from Normalization2.xhtml
    public String dataSummaryCheck() {
        RDataUtils.plotDataBox(sb, sb.getRConnection(), sb.getNewImage("qc_boxplot"), 72, "png");
        RDataUtils.plotDataBox2(sb, sb.getRConnection(), sb.getNewImage("qc_boxplot2"), 72, "png");
        
        RDataUtils.plotDataPCA(sb, sb.getRConnection(), sb.getNewImage("qc_pca"), 72, "png");
        RDataUtils.plotDataPCA2(sb, sb.getRConnection(), sb.getNewImage("qc_pca2"), 72, "png");
        
        RDataUtils.plotDataDensity(sb, sb.getRConnection(), sb.getNewImage("qc_density"), 72, "png");
        RDataUtils.plotDataDensity2(sb, sb.getRConnection(), sb.getNewImage("qc_density2"), 72, "png");
        
        RDataUtils.plotDataMeanStd(sb, sb.getRConnection(), sb.getNewImage("qc_meanstd"), 72, "png");
        RDataUtils.plotDataMeanStd2(sb, sb.getRConnection(), sb.getNewImage("qc_meanstd2"), 72, "png");
        
        RDataUtils.plotNormCheck(sb, sb.getRConnection(), sb.getNewImage("qc_normcheck"), 72, "png");
        RDataUtils.plotNormCheck2(sb, sb.getRConnection(), sb.getNewImage("qc_normcheck2"), 72, "png");
        sb.setImgDisabled(false);

        return "Data Summary";
    }    
    
    
    //Analysis start point!!! 
    public String enterAnalysis() {
        //first, need to reset those analyzed results for 
        String id = sb.getModuleType();
        if (id.equals("omics")) {
            return "OmicsAnalysis";
        } else {
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Failed to perform normalization!"));
            return null;
        }
    }
    
    //O2PLS Variable check command plots
    public void o2plsvariablecheck(int mode) {
        if (mode == 0){
            return;
        }
        RConnection RC = sb.getRConnection();
        String varOpt = analBean.getVarcriteriaOpt();
        
        if (varOpt.equals("adj")){ //cv adjusted R^2
            int res1 = RDataUtils.plotDatatable(sb, sb.getRConnection(), sb.getNewImage("adj_table"), 72, "png", jointc, orth1c, orth2c);
            int res2 = RDataUtils.plotDatatable2(sb, sb.getRConnection(), sb.getNewImage("adj_tableplot1"), 72, "png");
            int res3 = RDataUtils.plotDatatable3(sb, sb.getRConnection(), sb.getNewImage("adj_tableplot2"), 72, "png");
            int res4 = RDataUtils.plotDatatable4(sb, sb.getRConnection(), sb.getNewImage("adj_tableplot3"), 72, "png");
            
            if (res1 == 0 | res2 == 0 | res3 == 0 | res4 == 0) {
                String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
                FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            }
            
        } else if(varOpt.equals("pred")){//cv prediction mse
            int res1 = RDataUtils.plotDatatableprediction(sb, sb.getRConnection(), sb.getNewImage("prediction_table"), 72, "png", jointc, orth1c, orth2c);
            if (res1 == 0) {
                String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
                FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            }
        } else if(varOpt.equals("fixed")){ //fixed orth and joint
            int res1 = RDataUtils.plotDatatablefix1(sb, sb.getRConnection(), sb.getNewImage("jo2pre_orthfix"), 72, "png", joint, orth1, orth2, loops);
            int res2 = RDataUtils.plotDatatablefix2(sb, sb.getRConnection(), sb.getNewImage("jo1pre_orthfix"), 72, "png", joint, orth1, orth2, loops);
            int res3 = RDataUtils.plotDatatablefix3(sb, sb.getRConnection(), sb.getNewImage("orth2_joorth1fix"), 72, "png", joint, orth1, orth2, loops);
            int res4 = RDataUtils.plotDatatablefix4(sb, sb.getRConnection(), sb.getNewImage("orth1_joorth2fix"), 72, "png", joint, orth1, orth2, loops);
            
            if (res1 == 0 || res2 == 0 || res3 == 0 || res4 == 0) {
                String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
                FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            }            
            
            
        } else if(varOpt.equals("smpl")){ //pick # of orth and joint -> use joint, orth1, orth2
            int res1 = RDataUtils.plotDatatableuserpick(sb, sb.getRConnection(), sb.getNewImage("userpick_table"), 72, "png", joint, orth1, orth2);
            if (res1 == 0) {
                String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
                FacesContext.getCurrentInstance().addMessage(null,
                        new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
            }
            
        } else{
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", "Wrong option picked."));
        } 
    }
    
    
    //o2PLS main plots and tables
    public void o2plsrun() {
        //RConnection RC = sb.getRConnection();
        
        //First tab
        int res1 = RDataUtils.plotMaintable1(sb, sb.getRConnection(), sb.getNewImage("main_table1"), 72, "png", jointb, orth1b, orth2b);
        int res2 = RDataUtils.plotMaintable2(sb, sb.getRConnection(), sb.getNewImage("main_table2"), 72, "png", jointb, orth1b, orth2b);
        
        //Second tab
        int res3 = RDataUtils.plotMainplot1(sb, sb.getRConnection(), sb.getNewImage("main_plot1"), 72, "png", orth1b);
        int res4 = RDataUtils.plotMainplot2(sb, sb.getRConnection(), sb.getNewImage("main_plot2"), 72, "png", orth1b);
        
        //Third tab
        int res5 = RDataUtils.plotMainplot3(sb, sb.getRConnection(), sb.getNewImage("main_plot3"), 72, "png");   
        
        //Fourth tab
        int res6 = RDataUtils.plotMainplot4(sb, sb.getRConnection(), sb.getNewImage("main_plot4"), 72, "png", jointb); 
        
        //Fifth tab
        int res7 = RDataUtils.plotMainplot5(sb, sb.getRConnection(), sb.getNewImage("main_plot5"), 72, "png"); 
        int res8 = RDataUtils.plotMainplot6(sb, sb.getRConnection(), sb.getNewImage("main_plot6"), 72, "png"); 
        String res = RDataUtils.performQualityStat(sb, sb.getRConnection());
        setStats4prediction(res);
        
        //Sixth tab
        int res9 = RDataUtils.plotMainplot7(sb, sb.getRConnection(), sb.getNewImage("main_plot7"), 72, "png", jointb); 
        
        //Seventh tab
        int res10 = RDataUtils.plotMainplot8(sb, sb.getRConnection(), sb.getNewImage("main_plot8"), 72, "png"); 
        int res11 = RDataUtils.plotMainplot9(sb, sb.getRConnection(), sb.getNewImage("main_plot9"), 72, "png"); 
        int res12 = RDataUtils.plotMainplot10(sb, sb.getRConnection(), sb.getNewImage("main_plot10"), 72, "png");
        int res13 = RDataUtils.plotMainplot11(sb, sb.getRConnection(), sb.getNewImage("main_plot11"), 72, "png");
        
        if (res1 == 0 || res2 == 0 || res3 == 0 || res4 == 0 || res5 == 0 || res6 == 0 || res7 == 0 || res8 == 0 || res9 == 0 || res10 == 0 || res11 == 0 || res12 == 0 || res13 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }       
    } 
    
    //o2pls - plot SSQ update
    public void updateSSQTable() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.getSSQTable(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getTableChoice(), sb.getNewImage("main_tablessq"), 72, "png", orth1b, orth2b);
        int res2 = RDataUtils.plotMainplotupdate1(sb, RC, sb.getOmicsChoice(), sb.getNewImage("main_plot1"), 72, "png", orth1b, orth2b, sb.getTableN(), sb.getTableChoice());
        int res3 = RDataUtils.plotMainplotupdate2(sb, RC, sb.getOmicsChoice(), sb.getNewImage("main_plot2"), 72, "png", orth1b, orth2b, sb.getTableN(), sb.getTableChoice());
        
        if (res1 == 0 || res2 == 0 || res3 == 0 ) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    }    
    

    //o2pls - plot joint loading v.s. variable update
    public void updateLoadingTable() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.getLoadingTable(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getTableChoice(), sb.getNewImage("main_tableloading"), 72, "png", sb.getLoadingChoice(), sb.getLoadingChoice2(), sb.getLoadingChoice3(),jointb, orth1b, orth2b);
        int res2 = RDataUtils.plotMainplotupdate3(sb, RC, sb.getNewImage("main_plot3"), 72, "png", sb.getOmicsChoice(), sb.getLoadingChoice(), sb.getLoadingChoice2(), sb.getLoadingChoice3(), jointb, orth1b, orth2b, sb.getTableN(), sb.getTableChoice());
        int res3 = RDataUtils.getCorrelation(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getTableChoice(), sb.getNewImage("main_tableloadingcor"), 72, "png", sb.getLoadingChoice(), sb.getLoadingChoice2(), sb.getLoadingChoice3(),jointb, orth1b, orth2b, sb.getCorChoice());
        
        if (res1 == 0 || res2 == 0 || res3 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    }  
    
    
    //o2pls - plot joint loading v.s. loading update
    public void updateCompLoadingTable() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.getCompLoadingTable(sb, RC, sb.getOmicsChoice(), sb.getTableN(),sb.getNewImage("main_tableloadingcomp"), 72, "png", sb.getLoadingX(), sb.getLoadingY(), jointb);
        int res2 = RDataUtils.plotMainplotupdate4(sb, RC, sb.getNewImage("main_plot4"), 72, "png", sb.getOmicsChoice(), sb.getLoadingX(), sb.getLoadingY(), jointb, sb.getTableN());

        if (res1 == 0 || res2 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    }  
        
    //o2pls - quality prediction statistics string adj R^2
    private String stats4prediction = "";
    
    public String getStats4prediction() {
        return stats4prediction;
    }

    public void setStats4prediction(String stats4prediction) {
        this.stats4prediction = stats4prediction;
    }
    
    //final summary table 
    private String stats4prediction2 = "";
    
    public String getStats4prediction2() {
        return stats4prediction2;
    }

    public void setStats4prediction2(String stats4prediction2) {
        this.stats4prediction2 = stats4prediction2;
    }
    
    
    //PA stats
    private String stats4prediction3 = "";
    
    public String getStats4prediction3() {
        return stats4prediction3;
    }

    public void setStats4prediction3(String stats4prediction3) {
        this.stats4prediction3 = stats4prediction3;
    }
    
    //coinertia stat
    private String stats4prediction4 = "";
    
    public String getStats4prediction4() {
        return stats4prediction4;
    }

    public void setStats4prediction4(String stats4prediction4) {
        this.stats4prediction4 = stats4prediction4;
    }
    
    
    //PA update featrue contributions stats
    private String stats4prediction5 = "";
    
    public String getStats4prediction5() {
        return stats4prediction5;
    }

    public void setStats4prediction5(String stats4prediction5) {
        this.stats4prediction5 = stats4prediction5;
    }
    
    //coinertia update feature contributions stat
    private String stats4prediction6 = "";
    
    public String getStats4prediction6() {
        return stats4prediction6;
    }

    public void setStats4prediction6(String stats4prediction6) {
        this.stats4prediction6 = stats4prediction6;
    }
    
    
    
    
    
    //o2pls - plot quality of prediction least square update
    public void updatePredictionTable() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.getPredictionTable(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getJointPred(), sb.getNewImage("main_tableprediction"), 72, "png", jointb);
        int res2 = RDataUtils.plotMainplotupdate5(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getJointPred(), sb.getNewImage("main_plot5"), 72, "png", jointb, sb.getSampleMeta());
        int res3 = RDataUtils.plotMainplotupdate6(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getJointPred(), sb.getNewImage("main_plot6"), 72, "png", jointb);
        
        String res = RDataUtils.performQualityStatUpdate(sb, RC, sb.getOmicsChoice(), sb.getTableN(), sb.getJointPred(), jointb);
        setStats4prediction(res);
        
        if (res1 == 0 || res2 == 0 || res3 == 0 ) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    }            
    
    
    //o2pls - plot correlation of original and joint scores
    public void updateCorVarTable() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.getCorVarTable(sb, RC, sb.getOmicsChoice(), sb.getTableN(),sb.getNewImage("main_tablecorvar"), 72, "png", sb.getXxAxis(), sb.getYyAxis(), sb.getValueChoice(), jointb);
        int res2 = RDataUtils.plotMainplotupdate7(sb, RC, sb.getOmicsChoice(), sb.getTableN(),sb.getNewImage("main_plot7"), 72, "png", sb.getXxAxis(), sb.getYyAxis(), sb.getValueChoice(), sb.getVisualChoice(), jointb);
        
        if (res1 == 0 || res2 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    }      

    //o2pls - update correlation of joint scores and metadata
    public void updateCorrelation() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.plotMainplotupdate8(sb, RC, sb.getNewImage("main_plot8"), 72, "png", sb.getColorChoice(), sb.getShapeChoice());
        int res2 = RDataUtils.plotMainplotupdate9(sb, RC, sb.getNewImage("main_plot9"), 72, "png", sb.getColorChoice(), sb.getShapeChoice());
        int res3 = RDataUtils.plotMainplotupdate10(sb, RC, sb.getNewImage("main_plot10"), 72, "png", sb.getColorChoice(), sb.getShapeChoice());
        int res4 = RDataUtils.plotMainplotupdate11(sb, RC, sb.getNewImage("main_plot11"), 72, "png", sb.getColorChoice(), sb.getShapeChoice());
        
        if (res1 == 0 || res2==0 || res3==0 || res4==0) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    }  
    
    
    //o2pls - update summary table
    public void updateSummary() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.plotMainplotsummary(sb, RC, sb.getNewImage("main_plot12"), 72, "png", sb.getSampleChoice(), sb.getSampleNumber(), sb.getJointSummary(), jointb);
        int res2 = RDataUtils.plotMainplotsummary2(sb, RC, sb.getNewImage("main_plot13"), 72, "png", sb.getTopN(), sb.getTopN2(), sb.getJointSummary(), jointb);
        int res3 = RDataUtils.plotMainplotsummary2cor(sb, RC, sb.getNewImage("main_plot14crosscor"), 72, "png", sb.getTopN(), sb.getTopN2(), sb.getJointSummary(), jointb, sb.getCorChoiceb());
        
        String res = RDataUtils.performQualityStatUpdate2(sb, RC, sb.getSampleNumber(), sb.getJointSummary(), jointb, sb.getSampleChoice());
        setStats4prediction2(res);
        
        if (res1 == 0 | res2 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(RC);
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }
    } 
    
    
    
    
    
    
    //PA main plots and tables
    public void paRun() {
        RConnection RC = sb.getRConnection();

        int res1 = RDataUtils.plotPAplot1(sb, RC, sb.getNewImage("pa_plot1"), 72, "png", size, type, color, arrow); 
        int res2 = RDataUtils.plotPAplot2(sb, RC, sb.getNewImage("pa_plot2"), 72, "png", color); 
       
        String res = RDataUtils.performQualityStatUpdate3(sb, RC, perm);
        setStats4prediction3(res);
        
        String res0 = RDataUtils.performQualityStatUpdate4(sb, RC);
        setStats4prediction4(res0);
        
        if (res1 == 0 || res2 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }       
    } 
    
    
    //PA update feature contribution removing
    public void paRun2() {
        RConnection RC = sb.getRConnection();
        
        int res1 = RDataUtils.plotPAplotupdate1(sb, RC, sb.getNewImage("pa_plotupdate1"), 72, "png", choice, sb.getTopN1(), sb.getTopN2b(), adjustment, sb.getJoint(), sb.getOrth1(), sb.getOrth2(), sb.getJointChoice()); 
        int res2 = RDataUtils.plotPAplotupdate2(sb, RC, sb.getNewImage("pa_plotupdate2"), 72, "png", choice, sb.getTopN1(), sb.getTopN2b(), adjustment, sb.getJoint(), sb.getOrth1(), sb.getOrth2(), sb.getJointChoice()); 
        int res3 = RDataUtils.plotPAplotupdate3(sb, RC, sb.getNewImage("pa_plotupdate3"), 72, "png", choice, sb.getTopN1(), sb.getTopN2b(), adjustment, sb.getJoint(), sb.getOrth1(), sb.getOrth2(), sb.getJointChoice()); 

        String res4 = RDataUtils.performQualityStatUpdate5(sb, RC, choice, sb.getTopN1(), sb.getTopN2b(), adjustment, sb.getJoint(), sb.getOrth1(), sb.getOrth2(), sb.getJointChoice(), perm2);
        setStats4prediction5(res4);
        
        String res5 = RDataUtils.performQualityStatUpdate6(sb, RC, choice, sb.getTopN1(), sb.getTopN2b(), adjustment, sb.getJoint(), sb.getOrth1(), sb.getOrth2(), sb.getJointChoice());
        setStats4prediction6(res5);
        

        
        if (res1 == 0 || res2 == 0 || res3 == 0) {
            String errMsg = RDataUtils.getCurrentMsg(sb.getRConnection());
            FacesContext.getCurrentInstance().addMessage(null,
                    new FacesMessage(FacesMessage.SEVERITY_ERROR, "Error", errMsg));
        }       
    } 
    
    
    
    
}