package main;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import javax.sound.midi.Sequence;

import ngsep.sequences.DNASequence;
import ngsep.sequences.DNAShortKmer;
import ngsep.sequences.RawRead;
import ngsep.sequences.io.FastqFileReader;

public class ClusterTableManager {
	private static final char[] ALPH = {'A', 'C', 'G', 'T'};
	private static final String LOG = "/" + (new SimpleDateFormat("dd_HHmmss").format(new Date())) + "_" + "log.txt";
	private static String outFolder = "";
	private static int KMER_LENGTH = 31;
	private static int PREFIX = 10;
	private static String FILE_TYPE = ".fastq.gz";
	private static final boolean TEST_RUN = false;
	private static final int TABLE_SIZE = 2000000;
	private static final SequenceComparator SC = new SequenceComparator();
	
	private int newIndex = 0;
	private int[][] table;
	private HashMap<DNAShortKmer, Integer> index;
	private HashMap<DNAShortKmer, ArrayList<DNAShortKmer>> cluster;
	
	public ClusterTableManager(int size) {
		table = new int[4][size*KMER_LENGTH];
		index = new HashMap<>(size/2);
		if(TEST_RUN) cluster = new HashMap<>(size/2);
	}
	
	public void add(DNAShortKmer kmer) {
		boolean hit = false;
		for(int i = 0; i < kmer.length() && !hit; i++) {
			
			for(int j = 0; j < ALPH.length && !hit; j++) {
				Integer k = index.get(new DNAShortKmer(kmer.subSequence(0, i) + "" + ALPH[j] + kmer.subSequence(i+1, kmer.length())));
				if(k != null) {
					hit = true;
					append(kmer, k);
				}
			}
		}
		
		if(!hit) {
			newCluster(kmer);
		}
	}
	
	public void newCluster(DNAShortKmer kmer) {
		for(int i = 0; i < kmer.length(); i++) {
			char bp = Character.toUpperCase(kmer.charAt(i));
			if (!DNASequence.isInAlphabeth(bp)) {
				throw new RuntimeException("Non DNA kmer found " + kmer);
			}
			int ibp = DNASequence.BASES_STRING.indexOf(bp);
			table[ibp][newIndex*KMER_LENGTH + i]++;
		}
		
		index.put(kmer, newIndex++);
		if(TEST_RUN) {
			ArrayList<DNAShortKmer> clt = new ArrayList<>();
			clt.add(kmer);
			cluster.put(kmer, clt);
		}
	}
	
	public DNAShortKmer getRepresentative(int k) {
		char[] consensus = new char[KMER_LENGTH];
		
		for(int i = 0; i < consensus.length; i++) {
			int max = 0;
			for(int j = 0; j < DNASequence.BASES_STRING.length(); j++) {
				int next = table[j][k*KMER_LENGTH + i];
				if(max <= next) {
					consensus[i] = ALPH[j];
					max = next;
				}
					
			}
		}
		return new DNAShortKmer(new String(consensus));
	}
	
	public void append(DNAShortKmer kmer, int k) {
		DNAShortKmer oldKmer = getRepresentative(k);
		
		for(int i = 0; i < kmer.length(); i++) {
			boolean set = false;
			for(int j = 0; j < ALPH.length && !set; j++) {
				if(kmer.charAt(i) == ALPH[j]) {
					table[j][k*KMER_LENGTH + i]++;
					set = true;
				}
			}
			if(!set) System.err.println("ERROR AT " + kmer);
		}
		
		DNAShortKmer newKmer = getRepresentative(k);
		//System.out.println(oldKmer + " == " + newKmer);
		index.remove(oldKmer);
		index.put(newKmer, k);
		if(TEST_RUN) {
			ArrayList<DNAShortKmer> clt  = cluster.get(oldKmer);
			//System.out.println(oldKmer + " == " + clt);
			cluster.remove(oldKmer);
			clt.add(kmer);
			cluster.put(newKmer, clt);
		}
	}
	
	public void logClusterCount() {
		Set<DNAShortKmer> keys = index.keySet();
		
		try (PrintWriter log = new PrintWriter(new FileWriter(new File(outFolder + LOG), true))) {
			log.println("== START CLUSTER COUNT ==\n\n");
			for(DNAShortKmer kmer : keys) {
				int k = index.get(kmer);
				int count = 0;
				for(int i = 0; i < ALPH.length; i++) {
					count += table[i][k*KMER_LENGTH];
				}
				
				log.println(kmer + " == " + count);
			}
			log.println("\n\n== END CLUSTER COUNT ==\n\n");
			log.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void logClustersMembers() {
		if(!TEST_RUN) {
			System.err.println("Clusters members can only be logged on a test run."); return;
		}
		
		try (PrintWriter log = new PrintWriter(new FileWriter(new File(outFolder + LOG), true))) {
			log.println("== START CLUSTER MEMBERS ==\n\n");
			
			Set<DNAShortKmer> keys = cluster.keySet();
			for(DNAShortKmer kmer : keys) {
				log.println(kmer + "#");
				ArrayList<DNAShortKmer> members = cluster.get(kmer);
				for(DNAShortKmer member : members) log.println(member);
				log.println();
			}
			
			log.println("\n\n== END CLUSTER MEMBERS ==\n\n");
			log.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void process(String fileName) {
		int count = 0;
		try (FastqFileReader openFile = new FastqFileReader(fileName);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				//System.out.println(s);
				if(!s.contains("N")) {
					add(new DNAShortKmer(s.substring(PREFIX,PREFIX + KMER_LENGTH)));
					if((++count)%20000 == 0) {
						System.out.println("Processed " + (count) + " kmers");
						consistenceCheck();
					}
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("Processed a total of " + count + " kmers");
	}
	
	public void consistenceCheck() {
		for(int i = 0; i*KMER_LENGTH < table[0].length ; i++) {
			int count = 0;
			for(int k = 0; k < ALPH.length; k++) count += table[k][i*KMER_LENGTH];
			
			for(int j = 1; j < KMER_LENGTH; j++) {
				int pCount = 0;
				for(int k = 0; k < ALPH.length; k++) pCount += table[k][i*KMER_LENGTH + j];
				if(pCount != count) System.err.println("Count not adding up!");
			}
		}
		System.out.println("OKAY");
	}
	
	public void logDistancesInsideClustersDistr() {
		if(!TEST_RUN) {
			System.err.println("Cluster distances from the inside can only be logged on a test run."); return;
		}
		
		int[] freq = new int[KMER_LENGTH];
		
		try (PrintWriter log = new PrintWriter(new FileWriter(new File(outFolder + LOG), true))) {
			log.println("== START CLUSTER INSIDE DISTANCES ==\n\n");
			
			Set<DNAShortKmer> keys = cluster.keySet();
			for(DNAShortKmer kmer : keys) {
				int distance = 0;
				int count = 0;
				ArrayList<DNAShortKmer> members = cluster.get(kmer);
				for(int i = 0; i < members.size(); i++) {
					DNAShortKmer s1 = members.get(i);
					for(int j = i + 1; j < members.size(); j++) {
						DNAShortKmer s2 = members.get(j);
						distance += SC.compare(s1, s2);
						count++;
					}
				}
				
				try {
					int remainder = distance/count;
					freq[remainder]++;
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			
			log.println("Distance,Frequence");
			for(int i = 0; i < freq.length; i++) {
				log.println(i + "," + freq[i]);
			}
			
			log.println("\n\n== END CLUSTER INSIDE DISTANCES ==\n\n");
			log.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public void logClusterCountDistr() {
		int[] freq = new int[100000];
		Set<DNAShortKmer> keys = index.keySet();
		
		try (PrintWriter log = new PrintWriter(new FileWriter(new File(outFolder + LOG), true))) {
			log.println("== START CLUSTER COUNT ==\n\n");
			for(DNAShortKmer kmer : keys) {
				int k = index.get(kmer);
				int count = 0;
				for(int i = 0; i < ALPH.length; i++) {
					count += table[i][k*KMER_LENGTH];
				}
				
				freq[count]++;
			}
			
			log.println("Distance,Frequence");
			for(int i = 0; i < freq.length; i++) {
				log.println(i + "," + freq[i]);
			}
			
			log.println("\n\n== END CLUSTER COUNT ==\n\n");
			log.flush();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		ClusterTableManager cluster = new ClusterTableManager(TABLE_SIZE);
		cluster.process(args[0]);
		cluster.process(args[2]);
		outFolder = args[1];
		try {
			Thread.sleep(30000);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		cluster.logClusterCount();
		//cluster.logClustersMembers();
		cluster.logClusterCountDistr();
		//cluster.logDistancesInsideClustersDistr();
	}
}
