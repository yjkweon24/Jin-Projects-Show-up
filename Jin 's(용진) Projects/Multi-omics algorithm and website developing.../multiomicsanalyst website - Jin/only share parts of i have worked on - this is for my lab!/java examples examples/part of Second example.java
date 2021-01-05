import general.SessionBean1;
import org.rosuda.REngine.Rserve.RConnection;


public class RDataUtils {

    public static void initRData(RConnection RC) {
        try {
            RC.voidEval("Init.Data()");
        } catch (Exception rse) {
            System.out.println(rse);
        }
    }

    
    //help prcMsg table for datasummaryview.xhtml
    public static int[] sanityCheckData(RConnection RC, String datatype) {
        try {
            String rCommand = "SanityCheckData(\"" + datatype + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return (RC.eval(rCommand).asIntegers());
            
        } catch (Exception e) {
            System.out.println(e);
        }
        return null;
    }
    
    //help prcMsg table for datasummaryview.xhtml -> after missing value imputed
    public static int[] sanityCheckData3(RConnection RC, String datatype) {
        try {
            String rCommand = "SanityCheckData3(\"" + datatype + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return (RC.eval(rCommand).asIntegers());
            
        } catch (Exception e) {
            System.out.println(e);
        }
        return null;
    }
    
    //Missing value check
    public static int missingCheckData(RConnection RC, String datatype) {
        try {
            String rCommand = "MissingCheckData(\"" + datatype + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return (RC.eval(rCommand).asInteger());
            
        } catch (Exception e) {
            System.out.println(e);
        }
        return -99;
    }    
    
    //This is when the data is already normalized and filtered well -> For example datasets
    public static int[] CheckData(RConnection RC, String datatype) {
        try {
            String rCommand = "CheckData(\"" + datatype + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return (RC.eval(rCommand).asIntegers());

        } catch (Exception e) {
            System.out.println(e);
        }
        return null;
    }
    
    
    
    //First omics boxplot
    public static int plotDataBox(SessionBean1 sb, RConnection RC, String boxName, int dpi, String type) {
        try {
            String rCommand = "PlotDataBox(\"" + boxName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_boxplot", rCommand); // important in saving plots...

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    
    
    //Second omics boxplot
    public static int plotDataBox2(SessionBean1 sb, RConnection RC, String boxName, int dpi, String type) {
        try {
            String rCommand = "PlotDataBox2(\"" + boxName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_boxplot2", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }       
    
    
    //First omics PCA
    public static int plotDataPCA(SessionBean1 sb, RConnection RC, String pcaName, int dpi, String type) {
        try {
            String rCommand = "PlotDataPCA(\"" + pcaName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_pca", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    

    //Second omics PCA
    public static int plotDataPCA2(SessionBean1 sb, RConnection RC, String pcaName, int dpi, String type) {
        try {
            String rCommand = "PlotDataPCA2(\"" + pcaName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_pca2", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }        
    
    
    //First omics density
    public static int plotDataDensity(SessionBean1 sb, RConnection RC, String densityName, int dpi, String type) {
        try {
            String rCommand = "PlotDataDensity(\"" + densityName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_density", rCommand);
  
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    
    
    //Second omics density
    public static int plotDataDensity2(SessionBean1 sb, RConnection RC, String densityName, int dpi, String type) {
        try {
            String rCommand = "PlotDataDensity2(\"" + densityName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_density2", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }        
    
    
    
    //First omics density
    public static int plotDataMeanStd(SessionBean1 sb, RConnection RC, String meanStdName,int dpi, String type) {
        try {
            String rCommand = "PlotDataMeanStd(\"" + meanStdName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_meanstd", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }

    //Second omics density
    public static int plotDataMeanStd2(SessionBean1 sb, RConnection RC, String meanStdName,int dpi, String type) {
        try {
            String rCommand = "PlotDataMeanStd2(\"" + meanStdName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_meanstd2", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }

    //First omics normcheck
    public static int plotNormCheck(SessionBean1 sb, RConnection RC, String normName,int dpi, String type) {
        try {
            String rCommand = "PlotDataNorm(\"" + normName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_normcheck", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }

    //Second omics normcheck
    public static int plotNormCheck2(SessionBean1 sb, RConnection RC, String normName,int dpi, String type) {
        try {
            String rCommand = "PlotDataNorm2(\"" + normName + "\",\"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("qc_normcheck2", rCommand);

            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    
    
    
    
    //o2pls page variable pick adj R^2 table
    public static int plotDatatable(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int jointc, int orth1c, int orth2c) {
        try {
            String rCommand = "PlotDatatable(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + jointc + "\", \"" + orth1c + "\", \"" + orth2c + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("adj_table", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    } 
    
    public static int plotDatatable2(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotDatatable2(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("adj_tableplot1", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      
    
    public static int plotDatatable3(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotDatatable3(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("adj_tableplot2", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }     
    
    public static int plotDatatable4(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotDatatable4(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("adj_tableplot3", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }     
    
    
    
    

    //o2pls page variable pick prediction MSE table
    public static int plotDatatableprediction(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int jointc, int orth1c, int orth2c) {
        try {
            String rCommand = "PlotDatatableprediction(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + jointc + "\", \"" + orth1c + "\", \"" + orth2c + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("prediction_table", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
    
    //o2pls page variable pick by user (option 4) table
    public static int plotDatatableuserpick(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2) {
        try {
            String rCommand = "PlotDatatableuserpick(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("userpick_table", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
    //o2pls page variable pick by fixed (option3) - joint2 prediction orth fixed
    public static int plotDatatablefix1(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2, int loops) {
        try {
            String rCommand = "PlotDatatablefix1(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + loops + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("jo2pre_orthfix", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      
    
    // - joint1 prediction orth fixed
    public static int plotDatatablefix2(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2, int loops) {
        try {
            String rCommand = "PlotDatatablefix2(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + loops + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("jo1pre_orthfix", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      
     
    // orth1 joint fixed
    public static int plotDatatablefix3(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2, int loops) {
        try {
            String rCommand = "PlotDatatablefix3(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + loops + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("orth2_joorth1fix", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }       
     
    // orth2 joint fixed
    public static int plotDatatablefix4(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2, int loops) {
        try {
            String rCommand = "PlotDatatablefix4(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + loops + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("orth1_joorth2fix", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    
    
    // orth2 run - table 1 
     public static int plotMaintable1(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2) {
        try {
            String rCommand = "PlotMaintable1(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_table1", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }        
    
    // orth2 run - table 2 
     public static int plotMaintable2(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint, int orth1, int orth2) {
        try {
            String rCommand = "PlotMaintable2(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_table2", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }     
    
    // orth2 run - plot 1
     public static int plotMainplot1(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int orth1) {
        try {
            String rCommand = "PlotMainplot1(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + orth1 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot1", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      
    
    // orth2 run - plot 2
     public static int plotMainplot2(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int orth1) {
        try {
            String rCommand = "PlotMainplot2(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + orth1 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot2", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }       

     
    public static int getSSQTable(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int tableChoice, String imgName, int dpi, String format, int orth1, int orth2) {
        try {
            String rCommand = "PlotSSQTable(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + tableChoice + "\", \"" + orth1 + "\", \"" + orth2 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_tablessq", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }     
     
     
    // orth2 run - plot 1 update
     public static int plotMainplotupdate1(SessionBean1 sb, RConnection RC, int omicsChoice, String tableName, int dpi, String type, int orth1, int orth2, int tableN, int tableChoice) {
        try {
            String rCommand = "PlotMainplotupdate1(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + orth1 + "\", \"" + orth2 + "\" , \"" + omicsChoice + "\", \"" + tableN + "\", \"" + tableChoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot1", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      
    
    // orth2 run - plot 2 update
     public static int plotMainplotupdate2(SessionBean1 sb, RConnection RC, int omicsChoice, String tableName, int dpi, String type, int orth1, int orth2, int tableN, int tableChoice) {
        try {
            String rCommand = "PlotMainplotupdate2(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + orth1 + "\", \"" + orth2 + "\" , \"" + omicsChoice + "\", \"" + tableN + "\", \"" + tableChoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot2", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }           
    
    
    // orth2 run - plot 3
     public static int plotMainplot3(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot3(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot3", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }     
     
    // orth2 run - plot 3 table 
    public static int getLoadingTable(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int tableChoice, String imgName, int dpi, String format, int joint, int orth1, int orth2, int jointlimit, int orth1limit, int orth2limit) {
        try {
            String rCommand = "PlotLoadingTable(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + tableChoice + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + joint + "\", \"" + jointlimit + "\", \"" + orth1limit + "\", \"" + orth2limit + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_tableloading", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }          
     
    
    // orth2 run - plot 3 update 
    public static int plotMainplotupdate3(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int omicsChoice, int joint, int orth1, int orth2, int jointlimit, int orth1limit, int orth2limit, int tablen, int tablechoice) {
        try {
            String rCommand = "PlotMainplotupdate3(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\" , \"" + omicsChoice + "\", \"" + jointlimit + "\", \"" + orth1limit + "\", \"" + orth2limit + "\", \"" + tablen + "\", \"" + tablechoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot3", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }     
    
    // orth2 run - plot 3 correlation update 
    public static int getCorrelation(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int tableChoice, String imgName, int dpi, String format, int joint, int orth1, int orth2, int jointlimit, int orth1limit, int orth2limit, int corchoice) {
        try {
            String rCommand = "PlotUnivCor(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + tableChoice + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + joint + "\", \"" + jointlimit + "\", \"" + orth1limit + "\", \"" + orth2limit + "\", \"" + corchoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_tableloadingcor", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }    

    
    // orth2 run - plot 4  
    public static int plotMainplot4(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int joint) {
        try {
            String rCommand = "PlotMainplot4(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + joint + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot4", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    // orth2 run - plot 4 table 
    public static int getCompLoadingTable(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, String imgName, int dpi, String format, int loadingx, int loadingy, int jointb) {
        try {
            String rCommand = "PlotCompLoadingTable(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + loadingx + "\", \"" + loadingy + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_tableloadingcomp", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }         
    
    // orth2 run - plot 4 update 
    public static int plotMainplotupdate4(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int omicsChoice, int loadingx, int loadingy, int jointb, int tableN) {
        try {
            String rCommand = "PlotMainplotupdate4(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + loadingx + "\", \"" + loadingy + "\", \"" + omicsChoice + "\", \"" + jointb + "\", \"" + tableN + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot4", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      

    // orth2 run - plot 5
    public static int plotMainplot5(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot5(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot5", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    // orth2 run - plot 6 (table)
    public static int plotMainplot6(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot6(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot6", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   

    //o2pls run - plot6 table string
    public static String performQualityStat(SessionBean1 sb, RConnection RC) {
        try {
            String rCommand = "PerformQualityStat();";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   
    
    // orth2 run - plot 5 and 6 table
    public static int getPredictionTable(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int jointPred, String imgName, int dpi, String format, int jointb) {
        try {
            String rCommand = "PlotPredictionTable(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + jointPred + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_tableprediction", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }   

    // orth2 run - plot 5 update        
    public static int plotMainplotupdate5(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int jointPred, String imgName, int dpi, String format, int jointb, int samplemeta) {
        try {
            String rCommand = "PlotMainplotupdate5(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + jointPred + "\", \"" + jointb + "\", \"" + samplemeta + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot5", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      

    // orth2 run - plot 6 update        
    public static int plotMainplotupdate6(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int jointPred, String imgName, int dpi, String format, int jointb) {
        try {
            String rCommand = "PlotMainplotupdate6(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + jointPred + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot6", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      

    //o2pls run - plot6 table string update
    public static String performQualityStatUpdate(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, int jointPred, int jointb) {
        try {
            String rCommand = "PerformQualityStatUpdate(\"" + omicsChoice + "\", \"" + tableN + "\", \"" + jointPred + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   

    // orth2 run - plot 7
    public static int plotMainplot7(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int jointb) {
        try {
            String rCommand = "PlotMainplot7(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot7", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    

    // orth2 run - plot7 table update
    public static int getCorVarTable(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, String imgName, int dpi, String format, int xaxisjoint, int yaxisjoint, int option, int jointb) {
        try {
            String rCommand = "PlotCorVarTable(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + xaxisjoint + "\", \"" + yaxisjoint + "\", \"" + option + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_tablecorvar", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }            
    
    
    // orth2 run - plot7 plot update
    public static int plotMainplotupdate7(SessionBean1 sb, RConnection RC, int omicsChoice, int tableN, String imgName, int dpi, String format, int xaxisjoint, int yaxisjoint, int option, int visualChoice, int jointb) {
        try {
            String rCommand = "PlotMainplotupdate7(\"" + imgName + "\", \"" + dpi + "\", \"" + format + "\", \"" + omicsChoice + "\", \"" + tableN + "\", \"" + xaxisjoint + "\", \"" + yaxisjoint + "\", \"" + option + "\", \"" + visualChoice + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot7", rCommand);
            return (RC.eval(rCommand).asInteger());
        } catch (Exception e) {
        }
        return 0;
    }         
    
    // orth2 run - plot 8
    public static int plotMainplot8(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot8(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot8", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    // orth2 run - plot 9
    public static int plotMainplot9(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot9(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot9", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }       
    
    // orth2 run - plot 10
    public static int plotMainplot10(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot10(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot10", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    
    // orth2 run - plot 11
    public static int plotMainplot11(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type) {
        try {
            String rCommand = "PlotMainplot11(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot11", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    // orth2 run - plot 8 update
    public static int plotMainplotupdate8(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int color, int shape) {
        try {
            String rCommand = "PlotMainplotupdate8(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + color + "\", \"" + shape + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot8", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    // orth2 run - plot 9 update
    public static int plotMainplotupdate9(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int color, int shape) {
        try {
            String rCommand = "PlotMainplotupdate9(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + color + "\", \"" + shape + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot9", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    
    // orth2 run - plot 10 update
    public static int plotMainplotupdate10(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int color, int shape) {
        try {
            String rCommand = "PlotMainplotupdate10(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + color + "\", \"" + shape + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot10", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
    // orth2 run - plot 11 update
    public static int plotMainplotupdate11(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int color, int shape) {
        try {
            String rCommand = "PlotMainplotupdate11(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + color + "\", \"" + shape + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot11", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  

    // orth2 run - plot 2 summary
    public static int plotMainplotsummary(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int choice, int samples, int jointcomp, int jointb) {
        try {
            String rCommand = "PlotMainplotsummary(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + choice + "\", \"" + samples + "\", \"" + jointcomp + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot12", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
    public static int plotMainplotsummary2(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int tableNomics1, int tableNomics2, int jointcomp, int jointb) {
        try {
            String rCommand = "PlotMainplotsummary2(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + tableNomics1 + "\", \"" + tableNomics2+ "\", \"" + jointcomp + "\", \"" + jointb + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot13", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
        
    public static int plotMainplotsummary2cor(SessionBean1 sb, RConnection RC, String tableName, int dpi, String type, int tableNomics1, int tableNomics2, int jointcomp, int jointb, int correlationoption) {
        try {
            String rCommand = "PlotMainplotsummary2cor(\"" + tableName + "\", \"" + dpi + "\", \"" + type + "\", \"" + tableNomics1 + "\", \"" + tableNomics2+ "\", \"" + jointcomp + "\", \"" + jointb + "\", \"" + correlationoption + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("main_plot14crosscor", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
 
    
    public static String performQualityStatUpdate2(SessionBean1 sb, RConnection RC, int samplenumber, int jointchoice, int jointb, int choice) {
        try {
            String rCommand = "PerformQualityStatUpdate2(\"" + samplenumber + "\", \"" + jointchoice + "\", \"" + jointb + "\", \"" + choice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   
    
    //////////////////////
    
    // PA plot 1
    public static int plotPAplot1(SessionBean1 sb, RConnection RC, String imgName, int dpi, String type, float size, int types, int color, int arrow) {
        try {
            String rCommand = "PlotPAplot1(\"" + imgName + "\", \"" + dpi + "\", \"" + type + "\", \"" + size + "\", \"" + types + "\", \"" + color + "\", \"" + arrow + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("pa_plot1", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    // PA plot 1
    public static int plotPAplot2(SessionBean1 sb, RConnection RC, String imgName, int dpi, String type, int color) {
        try {
            String rCommand = "PlotPAplot2(\"" + imgName + "\", \"" + dpi + "\", \"" + type + "\", \"" + color + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("pa_plot2", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    public static String performQualityStatUpdate3(SessionBean1 sb, RConnection RC, int perm) {
        try {
            String rCommand = "PerformQualityStatUpdate3(\"" + perm + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   
    
    
    //Coinertia
    
    public static String performQualityStatUpdate4(SessionBean1 sb, RConnection RC) {
        try {
            String rCommand = "PerformQualityStatUpdate4();";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   
    
    
    
    

    //PA update - feature contributions
    public static int plotPAplotupdate1(SessionBean1 sb, RConnection RC, String imgName, int dpi, String type, int choice, int top1, int top2, int adjustment, int joint, int orth1, int orth2, int jointchoice) {
        try {
            String rCommand = "FeatureContributionPAplot1(\"" + imgName + "\", \"" + dpi + "\", \"" + type + "\", \"" + choice + "\", \"" + top1 + "\", \"" + top2 + "\", \"" + adjustment + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + jointchoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("pa_plotupdate1", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    public static int plotPAplotupdate2(SessionBean1 sb, RConnection RC, String imgName, int dpi, String type, int choice, int top1, int top2, int adjustment, int joint, int orth1, int orth2, int jointchoice) {
        try {
            String rCommand = "FeatureContributionPAplot2(\"" + imgName + "\", \"" + dpi + "\", \"" + type + "\", \"" + choice + "\", \"" + top1 + "\", \"" + top2 + "\", \"" + adjustment + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + jointchoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("pa_plotupdate2", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }   
    
    public static int plotPAplotupdate3(SessionBean1 sb, RConnection RC, String imgName, int dpi, String type, int choice, int top1, int top2, int adjustment, int joint, int orth1, int orth2, int jointchoice) {
        try {
            String rCommand = "FeatureContributionPAplot3(\"" + imgName + "\", \"" + dpi + "\", \"" + type + "\", \"" + choice + "\", \"" + top1 + "\", \"" + top2 + "\", \"" + adjustment + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + jointchoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            sb.addGraphicsCMD("pa_plotupdate3", rCommand); // important in saving plots...
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }  
    
    public static String performQualityStatUpdate5(SessionBean1 sb, RConnection RC, int choice, int top1, int top2, int adjustment, int joint, int orth1, int orth2, int jointchoice, int perm) {
        try {
            String rCommand = "PerformQualityStatUpdate5(\"" + choice + "\", \"" + top1 + "\", \"" + top2 + "\", \"" + adjustment + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + jointchoice + "\", \"" + perm + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   
    
    public static String performQualityStatUpdate6(SessionBean1 sb, RConnection RC, int choice, int top1, int top2, int adjustment, int joint, int orth1, int orth2, int jointchoice) {
        try {
            String rCommand = "PerformQualityStatUpdate6(\"" + choice + "\", \"" + top1 + "\", \"" + top2 + "\", \"" + adjustment + "\", \"" + joint + "\", \"" + orth1 + "\", \"" + orth2 + "\", \"" + jointchoice + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asString();
        } catch (Exception e) {
        }
        return null;
    }   
    
    
   
    
    //Get msgs from R function!!!
    public static String getCurrentMsg(RConnection RC) {
        try {
            return RC.eval("current.msg").asString();
        } catch (Exception e) {
        }
        return null;
    }
    
    //Get errors from R function!!!
    public static String getErrorMsg(RConnection RC) {
        try {
            return RC.eval("msg.vec").asString();
        } catch (Exception e) {
        }
        return null;
    }
    
   

    //Missing value imputaion for omics 1 
    public static int applyMissing(RConnection RC, String missImp) {
        try {
            String rCommand = "ApplyMissing1(\"" + missImp + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    

    //Missing value imputaion for omics 2 
    public static int applyMissing2(RConnection RC, String missImp) {
        try {
            String rCommand = "ApplyMissing2(\"" + missImp + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }    

    //Missing value imputaion for omics 1 and 2 
    public static int applyMissing3(RConnection RC, String missImp, String missImp2) {
        try {
            String rCommand = "ApplyMissing3(\"" + missImp + "\", \"" + missImp2 + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }      
    
    
    
    //Data filtering for the first omics
    public static int applyAbundanceFilter(RConnection RC, String filterOpt, int minCount, double filtlvl) {
        try {
            String rCommand = "ApplyAbundanceFilter(\"" + filterOpt + "\", " + minCount + ", " + filtlvl + ");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }

    public static int applyVarianceFilter(RConnection RC, String filtopt, double perct) {
        try {
            String rCommand = "ApplyVarianceFilter(\"" + filtopt + "\", " + perct + ");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }

    
    
        //Data filtering for the second omics
    public static int applyAbundanceFilter2(RConnection RC, String filterOpt, int minCount, double filtlvl) {
        try {
            String rCommand = "ApplyAbundanceFilter2(\"" + filterOpt + "\", " + minCount + ", " + filtlvl + ");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }

    public static int applyVarianceFilter2(RConnection RC, String filtopt, double perct) {
        try {
            String rCommand = "ApplyVarianceFilter2(\"" + filtopt + "\", " + perct + ");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }
    
    
    
    
    //Normalization for the first omics
    public static int performNormalization(RConnection RC, String scaleOpt, String transformOpt, float lambda1, float lambda1b, float lambda2b) {
        try {
            String rCommand = "PerformNormalization(\"" + scaleOpt + "\", \"" + transformOpt + "\", \"" + lambda1 + "\", \"" + lambda1b + "\", \"" + lambda2b + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }
    
    //Normalization for the second omics
    public static int performNormalization2(RConnection RC, String scaleOpt, String transformOpt, float lambda1, float lambda1b, float lambda2b) {
        try {
            String rCommand = "PerformNormalization2(\"" + scaleOpt + "\", \"" + transformOpt + "\", \"" + lambda1 + "\", \"" + lambda1b + "\", \"" + lambda2b + "\");";
            RCenter.recordRCommand(RC, rCommand);
            return RC.eval(rCommand).asInteger();
        } catch (Exception e) {
        }
        return 0;
    }
    
    
    
    //help prcMsg2 for summaryview
    public static int[] summarydata(RConnection RC) {
        try {
            String rCommand = "SanityCheckData2();";
            RCenter.recordRCommand(RC, rCommand);
            return (RC.eval(rCommand).asIntegers());
            
        } catch (Exception e) {
            System.out.println(e);
        }
        return null;
    }
    
    //Help on stack downloader!!!
    public static void saveAllData(RConnection RC, String analtype) {
        try {
            String rCommand = "SaveData(\"" + analtype + "\")";
            RCenter.recordRCommand(RC, rCommand);
            RC.voidEval(rCommand);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
