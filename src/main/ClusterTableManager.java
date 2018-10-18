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
	private static final int TABLE_SIZE = 2500000;
	private static final SequenceComparator SC = new SequenceComparator();
	
	private int newIndex = 0;
	private int kmerCount = 0;
	private int[][] table;
	private HashMap<DNAShortKmer, Integer> index;
	private HashMap<DNAShortKmer, ArrayList<DNAShortKmer>> cluster;
	
	public ClusterTableManager(int size) {
		table = new int[4][size*KMER_LENGTH];
		index = new HashMap<>(size/2);
		if(TEST_RUN) cluster = new HashMap<>(size/2);
	}
	
	/**
	 * Agrega el kmer por parametro al primer cluster a distancia de Hamming 1 o menos y registra sus datos en el consenso.
	 * Si el kmer no se relaciona con ningun cluster, crea un nuevo cluster cuyo consenso se compone unicamente de este kmer.
	 * @param kmer el kmer a agregar
	 */
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
	
	/**
	 * Genera un nuevo cluster cuyo unico miembro es el kmer por parametro.
	 * Reserva el espacio para el consenso e ingresa los datos del kmer.
	 * @param kmer representante del nuevo cluster
	 */
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
	
	/**
	 * Obtiene el kmer representativo del cluster cuyo concenso esta en la posicion k de la tabla de consensos
	 * @param k posicion en la tabla de consensos del cluster a evaluar
	 * @return
	 */
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
	
	/**
	 * Agrega la informacion del kmer al consenso del cluster que se encuentra en la posicion k de la tabla
	 * @param kmer el kmer a agregar al consenso
	 * @param k posicion en la tabla de consensos del cluster
	 */
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
	
	/**
	 * Cuenta la cantidad de miembros por cluster e imprime todo al log del programa 
	 */
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
	
	/**
	 * Imprime todos los miembros del cluster junto con su representativo (identificado por el prefijo #) al log del programa.
	 * Requiere que el programa corra en modo TEST_RUN
	 */
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
	
	/**
	 * Procesa todas las lecturas de un archivo fastq. Para cada lectura llama el metodo {@link #add(DNAShortKmer)}
	 * @param f archivo fastq para procesar
	 */
	public void process(File f) {
		try (FastqFileReader openFile = new FastqFileReader(f);) {
			Iterator<RawRead> reader = openFile.iterator();
			while(reader.hasNext()) {
				RawRead read = reader.next();
				String s = read.getSequenceString();
				//System.out.println(s);
				if(!s.contains("N")) {
					add(new DNAShortKmer(s.substring(PREFIX,PREFIX + KMER_LENGTH)));
					kmerCount++;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		System.out.println("\nProcessed a total of " + kmerCount + " kmers");
		consistenceCheck();
	}
	
	/**
	 * Revisa que no exista un concenso inconsistente.
	 * Un consenso inconsistente se considera como un consenso donde no todas las posiciones tienen la misma cantidad
	 * de votos.
	 */
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
	
	public void processDir(String path) {
		File[] files = (new File(path)).listFiles();
		for(File f : files) {
			if(f.getName().endsWith(FILE_TYPE)) process(f);
			System.out.println("Done with " + f.getName());
		}
	}
	
	/**
	 * Imprime las frecuencias de las distancias promedio para cada uno de los clusters de la tabla al log del programa.
	 * Requiere que el programa corra en modo TEST_RUN
	 */
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
	
	/**
	 * Imprime la frecuencia de la cantidad de elementos dentro de cada cluster al log del programa.
	 */
	public void logClusterCountDistr() {
		int[] freq = new int[1000000];
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
		long t = System.currentTimeMillis();
		ClusterTableManager cluster = new ClusterTableManager(TABLE_SIZE);
		cluster.processDir(args[0]);
		System.out.println(System.currentTimeMillis() - t);
		outFolder = args[1];
		
		cluster.logClusterCount();
		//cluster.logClustersMembers();
		cluster.logClusterCountDistr();
		//cluster.logDistancesInsideClustersDistr();
	}
}
