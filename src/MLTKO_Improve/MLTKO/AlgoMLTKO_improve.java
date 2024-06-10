package MLTKO_Improve.MLTKO;

import MLTKO.AlgoMLTKO;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryType;
import java.util.*;


public class AlgoMLTKO_improve {

    short version = 0x0104;
    /**
     * algorithm version number
     **/
    boolean useEUCP = false;
    /**
     * EUCP flag
     **/
    long candidateCount = 0;
    /**
     * candidate high-utility itemsets counter
     **/
    double totalTime = 0;
    /**
     * the time the algorithm terminated
     **/
    int huiCount = 0;
    /**
     * the number of HUI generated
     **/
    int K = 0;
    /**
     * the k parameter
     **/
    double minutility = 0;
    /**
     * the internal min utility variable
     **/

    public boolean debugging = false;        // for debugging purpose
    public boolean topdown = true;            // top-down mining or not
    public double algoRuntime = 0;        // recorded running time
    public double algoMemUsage = 0;        // recorded memory usage
    public boolean benchmark = false;        // benchmark flag, if true, the algorithm will only report runtime and memory

    PriorityQueue<ItemsetTKO> kItemsets;
    /**
     * the top k rules found until now
     **/

    long transCount = 0;
    /**
     * WIP: for scalability evaluations
     **/
    Taxonomy taxonomy = null;                /** for describing the taxonomy of a dataset **/

    /**
     * We create a map to store the TWU of each item
     */
    Map<Integer, Integer> mapItemToLevel;            // Item -> level hashmap
    Map<Integer, Double> mapItemToUtility;            // Item -> Utility
    Map<Integer, Double> mapItemToGWU;                // Item -> TWU/GWU
    Map<Integer, List<Integer>> mapItemToAncestor;    // Real taxonomy hashmap
    Map<Integer, Map<Integer, Double>> mapFMAP;        // EUCS:  key:item   key:another item   value:twu
    List<Set<Integer>> listItemHUI;
    Map<Integer, List<Integer>> mapItemToChild;

    /**
     * this class represent an item and its utility in a transaction
     */
    class Pair {
        int item = 0;
        /**
         * an item
         */
        double utility = 0;
        public Pair(){

        }
        public Pair(int item, double utility){
            this.item = item;
            this.utility = utility;
        }/** the utility of the item */
    }

    /**
     * Constructor
     */
    public AlgoMLTKO_improve(boolean eucp) {
        this.useEUCP = eucp;
    }

    /**
     * Run the algorithm
     *
     * @param input  the input file path
     * @param output the output file path
     * @param K      the parameter k
     * @throws IOException if an error occur for reading/writing to file.
     */
    public void runAlgorithm(String input, String tax, String output, int K)
            throws IOException {

        // allocate memory for the EUCP structure if requested in constructor
        if (this.useEUCP)
            mapFMAP = new HashMap<Integer, Map<Integer, Double>>();
        listItemHUI = new ArrayList<>();
        mapItemToGWU = new HashMap<Integer, Double>();
        mapItemToUtility = new HashMap<Integer, Double>();
        mapItemToLevel = new HashMap<Integer, Integer>();
        mapItemToAncestor = new HashMap<Integer, List<Integer>>();
        int ignore = 0;
        taxonomy = new Taxonomy(tax);
        mapItemToChild = new HashMap<>();
        //build mapItemToAncestor

        for (Map.Entry<Integer, Integer> entry : taxonomy.mapItemToParent.entrySet()) {
            Integer childItem = entry.getKey();
            Integer parentItem = entry.getValue();
            if (mapItemToAncestor.get(childItem) == null) {
                mapItemToAncestor.put(childItem, new ArrayList<>());
            }
            if (mapItemToAncestor.get(parentItem) == null) {
                mapItemToAncestor.put(parentItem, new ArrayList<>());
            }
            if (mapItemToChild.get(parentItem) == null) {
                mapItemToChild.put(parentItem, new ArrayList<>());

            }
            mapItemToChild.get(parentItem).add(childItem);
            List<Integer> ancestors = new ArrayList<>();
            Integer currentItem = childItem;
            while (currentItem != null) {
                Integer ancestor = taxonomy.mapItemToParent.get(currentItem);
                if (ancestor != null) {
                    ancestors.add(ancestor);
                }
                currentItem = ancestor;
            }
            mapItemToAncestor.put(childItem, ancestors);
        }
        int maxLevel = 0;
        for (Map.Entry<Integer, List<Integer>> entry : mapItemToAncestor.entrySet()) {
            mapItemToLevel.put(entry.getKey(), entry.getValue().size() + 1);
            if (maxLevel < entry.getValue().size() + 1) {
                maxLevel = entry.getValue().size() + 1;
            }
        }
        for (int i = 0; i < maxLevel; i++) {
            Set<Integer> list = new HashSet<>();
            listItemHUI.add(list);
        }

        long startTimestamp = System.currentTimeMillis();

        this.candidateCount = 0;
        this.minutility = 1;
        this.K = K;
        this.kItemsets = new PriorityQueue<ItemsetTKO>();

        // first dataset scan to calculate the TWU of each item.
        if (!benchmark) System.out.println("- First dataset scan...");

        Dataset dataset = new Dataset(input, Integer.MAX_VALUE);    // should perform similar transaction merging here, too
        transCount = dataset.getTransactions().size();

        for (int tid = 0; tid < transCount; tid++) {
            Transaction transaction = dataset.getTransactions().get(tid);
            double tu = transaction.getUtility();
            int[] items = transaction.getItems();
            double[] utils = transaction.getUtilities();
            HashMap<Integer, Double> ParentInTransaction = new HashMap<>();
            List<Set<Integer>> listAllAncestors = new ArrayList<>();
            for (int i = 0; i < maxLevel - 1; i++) {
                listAllAncestors.add(new HashSet<>());
            }
            for (int i = 0; i < items.length; i++) {        // for each item, add the transaction utility to its TWU
                Integer item = items[i];
                Double util = utils[i];
                Double twu = mapItemToGWU.get(item);// get the current TWU of that item

                for (int j = i + 1; j < items.length; j++) {
                    Integer itemJ = items[j];
                    Double utilityJ = utils[j];
                    Integer newItemI;
                    Integer newItemJ;
                    Double newUtilityI;
                    Double newUtilityJ;
                    if (item < itemJ) {
                        newItemI = item;
                        newItemJ = itemJ;
                        newUtilityI = util;
                        newUtilityJ = utilityJ;
                    } else {
                        newItemI = itemJ;
                        newItemJ = item;
                        newUtilityI = utilityJ;
                        newUtilityJ = util;
                    }
                    Map<Integer, Double> mapJ = mapFMAP.get(newItemI);
                    if (mapJ == null) {
                        mapJ = new HashMap<>();
                        mapFMAP.put(newItemI, mapJ);
                    }
                    Double utilityIJ = mapJ.get(newItemJ);
                    if (utilityIJ == null) {
                        mapJ.put(newItemJ, newUtilityI + newUtilityJ);
                    } else {
                        mapJ.put(newItemJ, newUtilityI + utilityIJ + newUtilityJ);
                    }
                }
                // add the utility of the item in the current transaction to its twu
                twu = (twu == null) ? tu : twu + tu;
                // calculate utility for each item on the dataset
                Double utility = mapItemToUtility.get(item);
                utility = (utility == null) ? util : utility + util;
                mapItemToUtility.put(item, utility);
                mapItemToGWU.put(item, twu);

                List<Integer> ancestors = mapItemToAncestor.get(item);
                for (Integer ancestorItem : ancestors) {
                    int level = mapItemToLevel.get(ancestorItem);
                    listAllAncestors.get(level - 1).add(ancestorItem);
                    Double utilityParent = ParentInTransaction.get(ancestorItem);
                    if (utilityParent == null) {
                        ParentInTransaction.put(ancestorItem, util);
                    } else {
                        ParentInTransaction.put(ancestorItem, utilityParent + util);
                    }
                }
            }
            for (int level = 0; level < listAllAncestors.size(); level++) {
                List<Integer> setItemInLevel = new ArrayList<>();
                setItemInLevel.addAll(listAllAncestors.get(level));
                setItemInLevel.sort(Collections.reverseOrder());
                for (int i = 0; i < setItemInLevel.size(); i++) {
                    Integer item = setItemInLevel.get(i);
                    Double utilityInTran = ParentInTransaction.get(item);
                    Double utility = mapItemToUtility.get(item);
                    Double twu = mapItemToGWU.get(item);
                    if (twu == null) {
                        mapItemToGWU.put(item, tu);
                        mapItemToUtility.put(item, utilityInTran);
                    } else {
                        mapItemToGWU.put(item, twu + tu);
                        mapItemToUtility.put(item, utility + utilityInTran);
                    }
                    Map<Integer, Double> mapJ = mapFMAP.get(item);
                    if (mapJ == null) {
                        mapJ = new HashMap<>();
                        mapFMAP.put(item, mapJ);
                    }
                    for (int j = i + 1; j < setItemInLevel.size(); j++) {
                        Integer itemJ = setItemInLevel.get(j);
                        Double utilityJ = ParentInTransaction.get(itemJ);
                        Double utilityIJ = mapJ.get(itemJ);
                        if (utilityIJ == null) {
                            mapJ.put(itemJ, utilityInTran + utilityJ);
                        } else {
                            mapJ.put(itemJ, utilityInTran + utilityIJ + utilityJ);
                        }
                    }
                }
            }
        } // for tid
        List<Double> allUtility = new ArrayList<>();
        for (Map.Entry<Integer, Map<Integer, Double>> entryF : mapFMAP.entrySet()) {
            Map<Integer, Double> mapJ = entryF.getValue();
            for (Map.Entry<Integer, Double> entryChild : mapJ.entrySet()) {
                allUtility.add(entryChild.getValue());
            }
        }
        for (Integer i : mapItemToUtility.keySet()) {
            allUtility.add(mapItemToUtility.get(i));
        }
        Collections.sort(allUtility, Collections.reverseOrder());

        // ==========================================================
        // ==== This section is based on Prof. Bay Vo suggestion ====

        //System.out.println(topUtils);
        if (K >= allUtility.size())
            minutility = 0;
        else
            minutility = allUtility.get(K - 1);
        if (minutility == 0) minutility = 1;
        if (!benchmark) System.out.println("  --> initial minUtil = " + minutility);
        if (!benchmark && debugging) {
            System.out.print("- Promising: ");
            for (Integer i : mapItemToUtility.keySet()) {
                System.out.print(i + " ");
            }
            System.out.println();
        }

        // ==========================================================
        // for faster accessing the utility lists, they are are stored using map as pair:
        // <KEY: item, VALUE: utility list associated to that item>
        List<List<List<Pair>>> revisedTransactions = new ArrayList<>();
        for (int i = 0;i < maxLevel; i++) {

            List<List<Pair>> revisedTransactionTemp = new ArrayList<>();

            for (int j = 0; j < transCount; j++) {
                List<Pair> rrTemp = new ArrayList<Pair>();
                List<Integer> ctTemp = new ArrayList<Integer>();
                revisedTransactionTemp.add(rrTemp);
            } // for j

            revisedTransactions.add(revisedTransactionTemp);
        } // for i

        if (!this.benchmark) {
            System.out.println("==== DATASET CHARACTERISTICS ====");
            System.out.println(" Dataset: <" + input + " / " + tax + ">");
            System.out.println(" |D|    : " + transCount);
            System.out.println(" |GI|   : " + taxonomy.parentCount());
            System.out.println(" Depth  : " + maxLevel);
            System.out.println(" T_max  : " + dataset.getMaxTransLength());
            System.out.println(" T_avg  : " + dataset.getAvgTransLength());
            System.out.println("=================================");

            System.out.println("- Second dataset scan...");
        }
        mapFMAP = new HashMap<>();

        for (int tid = 0; tid < transCount; tid++) {
            Transaction transaction = dataset.getTransactions().get(tid);
            int[] items = transaction.getItems();
            double[] utils = transaction.getUtilities();
            Map<Integer,Double> parentInTransaction = new HashMap<>();
            for (int i = 0; i < items.length; i++) {
                Integer item = items[i];
                if(mapItemToGWU.get(item)>=minutility){
                    Pair pair = new Pair();
                    pair.item = item;
                    pair.utility = utils[i];
                    int level = mapItemToLevel.get(item);
                    revisedTransactions.get(level - 1).get(tid).add(pair);
                }
                List<Integer> listAllAncestors = mapItemToAncestor.get(item);
                for (Integer itemParent : listAllAncestors) {
                    if(mapItemToGWU.get(itemParent)>=minutility){
                        Double utilityParent = parentInTransaction.get(itemParent);
                        if(utilityParent==null){
                            parentInTransaction.put(itemParent,utils[i]);
                        }else{
                            parentInTransaction.put(itemParent,utilityParent+utils[i]);
                        }
                    }
                }
            }
            for(Map.Entry<Integer,Double> entry : parentInTransaction.entrySet()){
                int level = mapItemToLevel.get(entry.getKey())-1;
                revisedTransactions.get(level).get(tid).add(new Pair(entry.getKey(),entry.getValue()));
            }

            for (int i = 0; i < maxLevel; i++) {                // sort the transactions
                Collections.sort(revisedTransactions.get(i).get(tid), new Comparator<Pair>() {
                    public int compare(Pair o1, Pair o2) {
                        return compareItems(o1.item, o2.item);
                    }
                });
            }
        }

        // BUild UL
        mapFMAP = new HashMap<Integer, Map<Integer, Double>>();
        Map<Integer, UtilityList> mapItemToUtilityList = new HashMap<Integer, UtilityList>();

        for (int level = 0; level < maxLevel; level++) {
            System.out.println(kItemsets.size());
            System.out.println(minutility);
            if (!benchmark) System.out.println("- Constructing utility lists for level" + level);
            List<UtilityList> UtilityListOfILevel = new ArrayList<>();
            if (level == 0) {
                for (Integer item : mapItemToGWU.keySet()) {
                    if (mapItemToGWU.get(item) >= this.minutility && mapItemToLevel.get(item) == level + 1) {
                        UtilityList uList = new UtilityList(item);
                        mapItemToUtilityList.put(item, uList);
                        UtilityListOfILevel.add(uList);
                    }
                }
            } else {
                List<Integer> ChildOfPreLevel = new ArrayList<>();
                Set<Integer> ItemsOfLevel = listItemHUI.get(level - 1);
                for (Integer parent : ItemsOfLevel) {
                    List<Integer> listChild = mapItemToChild.get(parent);
                    if (listChild != null) {
                        for (Integer child : listChild) {
                            ChildOfPreLevel.add(child);
                        }
                    }
                }
                for (Integer item : ChildOfPreLevel) {
                    if (mapItemToGWU.get(item) != null && mapItemToGWU.get(item) >= this.minutility) {
                        UtilityList uList = new UtilityList(item);
                        mapItemToUtilityList.put(item, uList);
                        UtilityListOfILevel.add(uList);
                    }
                }
            }
            Collections.sort(UtilityListOfILevel, (o1, o2) -> {
                return compareItems(o1.item, o2.item);        // compare the TWU of the items
            });

            List<List<Pair>> DBinLevel = revisedTransactions.get(level);
            for (int tid = 0; tid < DBinLevel.size(); tid++) {
                List<Pair> transaction = DBinLevel.get(tid);
                double remainingUtility = 0;
                double newTWU = 0;
                List<Pair> newTran = new ArrayList<>();
                for (int i = transaction.size() - 1; i >= 0; i--) {
                    Pair pair = transaction.get(i);
                    UtilityList utilityListOfItem = mapItemToUtilityList.get(pair.item);
                    if (utilityListOfItem != null) {
                        Element element = new Element(tid, pair.utility, remainingUtility);
                        utilityListOfItem.addElement(element);
                        remainingUtility += pair.utility;
                        newTWU += pair.utility;
                        newTran.add(0, pair);
                    }
                }
                for (int i = 0; i < newTran.size(); i++) {
                    Pair pair = newTran.get(i);
                    Map<Integer, Double> mapFMAPItem = mapFMAP.get(pair.item);
                    if (mapFMAPItem == null) {
                        mapFMAPItem = new HashMap<Integer, Double>();
                        mapFMAP.put(pair.item, mapFMAPItem);
                    }
                    for (int j = i + 1; j < newTran.size(); j++) {
                        Pair pairAfter = newTran.get(j);
                        Double twuSum = mapFMAPItem.get(pairAfter.item);
                        if (twuSum == null) {
                            mapFMAPItem.put(pairAfter.item, newTWU);
                        } else {
                            mapFMAPItem.put(pairAfter.item, twuSum + newTWU);
                        }
                    }
                }
            }
            if (!benchmark) System.out.println("- MLHUI mining...");
            search(new int[0], null, UtilityListOfILevel);
        }

        if (!benchmark && useEUCP)
            if (debugging) {
                for (Integer i : mapItemToUtility.keySet()) {
                    HashMap<Integer, Double> mapFMAPItem = (HashMap<Integer, Double>) mapFMAP.get(i);
                    for (Integer j : mapItemToUtility.keySet()) {
                        Double v = mapFMAPItem.get(j);
                        if (v != null)
                            System.out.print("[" + i + "," + j + "]=" + v + " ");
                    }
                    System.out.println();
                }
            }

        // dataset and taxonomy is now no longer needed, discard them
        dataset = null;
        taxonomy = null;
        if (output != null) {
            System.out.println("- Saving results...");
            writeResultTofile(output);
        }


        totalTime = (System.currentTimeMillis() - startTimestamp) / 1000.0;

        algoRuntime = totalTime;
        algoMemUsage =

                peakHeapUsage();

        if (!benchmark) System.out.println("- Done.");
    }

    /**
     * This is the recursive method to find all high utility itemsets. It writes the itemsets to the output file.
     *
     * @param prefix This is the current prefix. Initially, it is empty.
     * @param pUL    This is the Utility List of the prefix. Initially, it is empty.
     * @param ULs    The utility lists corresponding to each extension of the prefix.
     * @throws IOException
     */
    private void search(int[] prefix, UtilityList pUL, List<UtilityList> ULs) throws IOException {

        // For each extension X of prefix P
        for (int i = 0; i < ULs.size(); i++) {
            UtilityList X = ULs.get(i);

            // If pX is a high utility itemset, we save the itemset: pX
            if (X.sumIutils >= minutility)
                writeOut(prefix, X.item, X.sumIutils);

            // Z-item pruning
            if (X.sumRutils == 0) continue;

            // If sum remaining utilities of pX >= minUtility, we explore extensions of pX.
            // Fig 5, Line 06 from the TKO paper (this is the pruning condition)
            if (X.sumRutils + X.sumIutils >= minutility) {

                // This list will contain the utility lists of pX extensions.
                List<UtilityList> exULs = new ArrayList<UtilityList>();

                // For each extension of p appearing after X according to the ascending order
                for (int j = i + 1; j < ULs.size(); j++) {
                    UtilityList Y = ULs.get(j);

                    // ==== EUCP ====
                    if (this.useEUCP) {
                        Map<Integer, Double> mapTWUF = mapFMAP.get(X.item);
                        if (mapTWUF != null) {
                            Double twuF = mapTWUF.get(Y.item);
                            if (twuF == null || twuF < minutility) continue;
                        }
                    }
                    // ==============

                    candidateCount++;

                    // we construct the extension pXY and add it to the list of extensions of pX
                    exULs.add(construct(pUL, X, Y));
                }

                // We create new prefix pX
                int[] newPrefix = new int[prefix.length + 1];
                System.arraycopy(prefix, 0, newPrefix, 0, prefix.length);
                newPrefix[prefix.length] = X.item;

                // We make a recursive call to discover all itemsets with the prefix pX
                search(newPrefix, X, exULs);
            }

        }
    }

    /**
     * Method to dynamically increase minutil and maintain the list of top-K MLHUIs
     *
     * @param prefix  itemset
     * @param item    to be appended to the prefix
     * @param utility the utility of the prefix concatenated with the item
     */
    private void writeOut(int[] prefix, int item, double utility) {
        ItemsetTKO itemset = new ItemsetTKO(prefix, item, utility);
        kItemsets.add(itemset);
        if (utility > this.minutility) {
            int level = mapItemToLevel.get(item) - 1;
            for (int itemInItemset : prefix) {
                listItemHUI.get(level).add(itemInItemset);
            }
            listItemHUI.get(level).add(item);
        }
        if (kItemsets.size() > this.K) {
            if (utility > this.minutility) {
                ItemsetTKO lower;
                do {
                    lower = kItemsets.peek();
                    if (lower == null) break; // / IMPORTANT
                    kItemsets.remove(lower);
                } while (kItemsets.size() > this.K);
                this.minutility = kItemsets.peek().utility;
                //System.out.println(String.format("%.5f", this.minutility));
            }
        }
    }

    /**
     * This method constructs the utility list of pXY
     *
     * @param P  :  the utility list of prefix P.
     * @param px : the utility list of pX
     * @param py : the utility list of pY
     * @return the utility list of pXY
     */
    private UtilityList construct(UtilityList P, UtilityList px, UtilityList py) {
        UtilityList pxyUL = new UtilityList(py.item);    // create an empy utility list for pXY
        int xSize = px.elements.size();

        for (int m = 0; m < xSize; m++) {    // for each element in the utility list of pX
            Element ex = px.elements.get(m);

            // find element ey in py with tid = ex.tid
            Element ey = findElementWithTID(py, ex.tid);
            if (ey == null) continue;

            if (P == null) {    // if the prefix p is null
                Element eXY = new Element(ex.tid,    // Create the new element
                        ex.iutils + ey.iutils,
                        ey.rutils);
                pxyUL.addElement(eXY);                    // add the new element to the utility list of pXY

            } else {
                // find the element in the utility list of p wih the same tid
                Element e = findElementWithTID(P, ex.tid);
                if (e != null) {
                    Element eXY = new Element(ex.tid,    // Create new element
                            ex.iutils + ey.iutils - e.iutils,
                            ey.rutils);
                    pxyUL.addElement(eXY);                // add the new element to the utility list of pXY
                }
            }
        }
        return pxyUL;                                    // return the utility list of pXY.
    }

    /**
     * Do a binary search to find the element with a given tid in a utility list
     *
     * @param ulist the utility list
     * @param tid   the tid
     * @return the element or null if none has the tid.
     */
    private Element findElementWithTID(UtilityList ulist, int tid) {
        List<Element> list = ulist.elements;

        int first = 0, last = list.size() - 1;
        while (first <= last) {
            int m = (first + last) >>> 1;
            if (list.get(m).tid < tid) first = m + 1;
            else if (list.get(m).tid > tid) last = m - 1;
            else
                return list.get(m);
        }
        return null;
    }

    /**
     * Write the result to a file
     *
     * @param path the output file path
     * @throws IOException if an exception for reading/writing to file
     */
    public void writeResultTofile(String path) throws IOException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(path));
        Iterator<ItemsetTKO> iter = kItemsets.iterator();
        while (iter.hasNext()) {
            StringBuffer buffer = new StringBuffer();
            ItemsetTKO itemset = (ItemsetTKO) iter.next();

            // append the prefix
            for (int i = 0; i < itemset.getItemset().length; i++) {
                buffer.append(itemset.getItemset()[i]);
                buffer.append(' ');
            }
            buffer.append(itemset.item);

            // append the utility value
            buffer.append(" #UTIL: ");
            buffer.append(String.format("%.5f", itemset.utility));

            // write to file
            writer.write(buffer.toString());
            if (iter.hasNext()) {
                writer.newLine();
            }
        }
        writer.close();
    }

    private int compareItems(int item1, int item2) {
        int compare = (int) (mapItemToGWU.get(item1) - mapItemToGWU.get(item2));
        // if the same, use the lexical order otherwise use the TWU
        return (compare == 0) ? item1 - item2 : compare;
    }

    public double peakHeapUsage()    // return the peak heap usage.
    {
        double retVal = 0;
        try {
            List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
            double total = 0;
            for (MemoryPoolMXBean memoryPoolMXBean : pools) {
                if (memoryPoolMXBean.getType() == MemoryType.HEAP) {
                    long peakUsed = memoryPoolMXBean.getPeakUsage().getUsed();
                    total = total + peakUsed;
                }
            }
            retVal = total / (1024L * 1024L);
        } catch (Throwable t) {
            System.err.println("Exception in agent: " + t);
        }
        return retVal;
    }

    /**
     * Print statistics about the latest execution to System.out.
     */
    public void printStats() {
        System.out.println("=============  MLTKO IMPROVE - " + getVersion() + " =============");
        System.out.println(" K           : " + this.K);
        System.out.println(" EUCP enabled: " + (this.useEUCP ? "yes" : "no"));
        System.out.println(" Top-down    : " + (this.topdown ? "yes" : "no"));
        System.out.println(" MLHUIs      : " + kItemsets.size());
        System.out.println(" Candidates  : " + candidateCount);
        System.out.println(" Runtime     : ~" + String.format("%.3f", totalTime) + " s");
        System.out.println(" Peak memory : ~" + String.format("%.3f", peakHeapUsage()) + " MB");
        System.out.println("==========================================");
    }

    /**
     * most useless stuffs evah
     **/
    public String getVersion() {
        return "v." + ((version & 0xFF00) >> 8) + "." + (version & 0xFF);
    }

    public void printLogo() {
        System.out.println("       _     _____ _____ _____");
        System.out.println(" _____| |___|_   _|  |  |     |");
        System.out.println("|     | |___| | | |    -|  |  |");
        System.out.println("|_|_|_|_|     |_| |__|__|_____|");
        System.out.println("                          " + getVersion());
    }
}
