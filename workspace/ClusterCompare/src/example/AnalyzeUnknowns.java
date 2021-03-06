package example;

import org.w3c.dom.*;
import javax.xml.parsers.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class AnalyzeUnknowns 
{

	public static void main(String[] args) 
	{
		int PPMTolerance = 10;
		// TODO Auto-generated method stub
		try{
			//Generates a List of each database compounds
			List<Element> HMDBList = new ArrayList<>();
			List<Element> KEGGList = new ArrayList<>();
			List<Element> NISTList = new ArrayList<>();
			List<Element> PUBCHEMList = new ArrayList<>();
			
			
			//For each database it has to load in the compounds and populate the lists from the XML file
			//that is associated with it
			File keggInputFile = new File("/Users/MatthewDeitz/Desktop/MetaboliteDB/KEGGClusterResults.xml");
			DocumentBuilderFactory keggdbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder keggdBuilder = keggdbFactory.newDocumentBuilder();
			Document keggdoc = keggdBuilder.parse(keggInputFile);
			keggdoc.getDocumentElement().normalize();
			NodeList keggnList = keggdoc.getElementsByTagName("row");
			for(int temp = 0; temp < keggnList.getLength();temp++)
			{
				Element e = ((Element) keggnList.item(temp));
				if (!e.getElementsByTagName("CenterMass").item(0).getTextContent().equals("NA"))
				{
					KEGGList.add(e);
				}
			}
			
			File hmdbInputFile = new File("/Users/MatthewDeitz/Desktop/MetaboliteDB/HMDBClusters.xml");
			DocumentBuilderFactory hmdbdbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder hmdbdBuilder = hmdbdbFactory.newDocumentBuilder();
			Document hmdbdoc = hmdbdBuilder.parse(hmdbInputFile);
			hmdbdoc.getDocumentElement().normalize();
			NodeList hmdbnList = hmdbdoc.getElementsByTagName("row");
			for(int temp = 0; temp < hmdbnList.getLength();temp++)
			{
				Element e = ((Element) hmdbnList.item(temp));
				if (!e.getElementsByTagName("CenterMass").item(0).getTextContent().equals("NA"))
				{
					HMDBList.add(e);
				}
			}
			
			File nistInputFile = new File("/Users/MatthewDeitz/Desktop/MetaboliteDB/NISTClusters.xml");
			DocumentBuilderFactory nistdbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder nistdBuilder = nistdbFactory.newDocumentBuilder();
			Document nistdoc = nistdBuilder.parse(nistInputFile);
			nistdoc.getDocumentElement().normalize();
			NodeList nistnList = nistdoc.getElementsByTagName("row");
			for(int temp = 0; temp < nistnList.getLength();temp++)
			{
				Element e = ((Element) nistnList.item(temp));
				if (!e.getElementsByTagName("CenterMass").item(0).getTextContent().equals("NA"))
				{
					NISTList.add(e);
				}
			}
			
			
			File pubInputFile = new File("/Users/MatthewDeitz/Desktop/MetaboliteDB/PUBCHEMlusterResults.xml");
			DocumentBuilderFactory pubdbFactory = DocumentBuilderFactory.newInstance();
			DocumentBuilder pubdBuilder = pubdbFactory.newDocumentBuilder();
			Document pubdoc = pubdBuilder.parse(pubInputFile);
			pubdoc.getDocumentElement().normalize();
			NodeList pubnList = pubdoc.getElementsByTagName("row");
			for(int temp = 0; temp < pubnList.getLength();temp++)
			{
				Element e = ((Element) pubnList.item(temp));
				if (!e.getElementsByTagName("CenterMass").item(0).getTextContent().equals("NA"))
				{
					PUBCHEMList.add(e);
				}
			}
			
			
			
			
			//read in the unknown compounds to be matched to the databases into a list
			List <String> unknowns = new ArrayList<>();
			BufferedReader br = new BufferedReader(new FileReader("/Users/MatthewDeitz/Desktop/test.txt"));
			String line = "";
			while ((line = br.readLine())!=null)
			{
				unknowns.add(line.trim());
			}
			br.close();
			
			
			//attempt to open the results file and paste if it has been found
			FileWriter write = new FileWriter("/Users/MatthewDeitz/Desktop/results.txt");
			PrintWriter results = new PrintWriter(write);
			for(String unknown : unknowns)//for all unknowns
			{
				double unknownmz = Double.parseDouble(unknown);//parse the mz value into a double
				boolean found = false;
				//format for searching every database
				for(int temp = 0;temp < HMDBList.size();temp++)
				{
					Element HMDBMetabolite = HMDBList.get(temp);//get the element in the list		
					Double HMDBmz = Double.parseDouble(HMDBMetabolite.getElementsByTagName("CenterMass").item(0).getTextContent().trim());//compare it with the center mass of the cluster
					if ((HMDBmz-unknownmz)/unknownmz*1000000 < -PPMTolerance)
					{
						//do nothing because not within tolerance
					}
					else if ((HMDBmz-unknownmz)/unknownmz*1000000 > PPMTolerance)
					{
						//do nothing because not within tolerance
					}
					else
					{
						//matches the database
						String exactMass = HMDBMetabolite.getElementsByTagName("Monisotopic_weight").item(0).getTextContent().trim();
						if(exactMass.indexOf(',')<0)
						{
							//if there is only one compound in the cluster
							results.println(unknown+" www.hmdb.ca/metabolites/"+HMDBMetabolite.getElementsByTagName("Accession").item(0).getTextContent());
						}
						else //if it is matched with a cluster then find closest match in that cluster
						{
							String[] parts = HMDBMetabolite.getElementsByTagName("Monisotopic_weight").item(0).getTextContent().substring(1).split(", ");
							int index = 0;
							int closestIndex = 0;
							Double closestVal = 1.0;
							for(String clustermz : parts)
							{
								Double mzValue = Double.parseDouble(clustermz);
								if (mzValue - unknownmz < closestVal)
									{
										closestVal = mzValue;
										closestIndex = index;
									}
								index++;
							}
							results.println(unknown+" www.hmdb.ca/metabolites/"+HMDBMetabolite.getElementsByTagName("Accession").item(0).getTextContent().split(", ")[closestIndex]);
							
						}

						found = true;
					}
					
				}
				
				for(int temp = 0; temp < KEGGList.size();temp++)
				{
					Element KEGGMetabolite = KEGGList.get(temp);
					Double KEGGmz = Double.parseDouble(KEGGMetabolite.getElementsByTagName("CenterMass").item(0).getTextContent().trim());
					if ((KEGGmz-unknownmz)/unknownmz*1000000 < -PPMTolerance)
					{
						//do nothing
					}
					else if ((KEGGmz-unknownmz)/unknownmz*1000000 > PPMTolerance)
					{
						//do nothing
					}
					else
					{
						String exactMass = KEGGMetabolite.getElementsByTagName("Exact_Mass").item(0).getTextContent().trim();
						if(exactMass.indexOf(',')<0)
						{
							results.println(unknown+" http://www.genome.jp/dbget-bin/www_bget?cpd:"+KEGGMetabolite.getElementsByTagName("Entry").item(0).getTextContent());
						}
						else
						{
							String[] parts = KEGGMetabolite.getElementsByTagName("Exact_Mass").item(0).getTextContent().substring(1).split(", ");
							int index = 0;
							int closestIndex = 0;
							Double closestVal = 1.0;
							for(String clustermz : parts)
							{
								Double mzValue = Double.parseDouble(clustermz);
								if (mzValue - unknownmz < closestVal)
									{
										closestVal = mzValue;
										closestIndex = index;
									}
								index++;
							}
							results.println(unknown+" http://www.genome.jp/dbget-bin/www_bget?cpd:"+KEGGMetabolite.getElementsByTagName("Entry").item(0).getTextContent().split(", ")[closestIndex]);
							
						}
						
						found = true;
					}
					
				}
				
				for(int temp = 0; temp < NISTList.size();temp++)
				{
					Element NISTMetabolite = NISTList.get(temp);
					Double NISTmz = Double.parseDouble(NISTMetabolite.getElementsByTagName("CenterMass").item(0).getTextContent().trim());
					if ((NISTmz-unknownmz)/unknownmz*1000000 < -PPMTolerance)
					{
						//do nothing
					}
					else if ((NISTmz-unknownmz)/unknownmz*1000000 > PPMTolerance)
					{
						//do nothing
					}
					else
					{
						String exactMass = NISTMetabolite.getElementsByTagName("ExactMass").item(0).getTextContent().trim();
						if(exactMass.indexOf(',')<0)
						{
							results.println(unknown+" https://www.ncbi.nlm.nih.gov/pccompound/?term="+NISTMetabolite.getElementsByTagName("InChIKey").item(0).getTextContent());
						}
						else
						{
							String[] parts = NISTMetabolite.getElementsByTagName("ExactMass").item(0).getTextContent().substring(1).split(", ");
							int index = 0;
							int closestIndex = 0;
							Double closestVal = 1.0;
							for(String clustermz : parts)
							{
								Double mzValue = Double.parseDouble(clustermz);
								if (mzValue - unknownmz < closestVal)
									{
										closestVal = mzValue;
										closestIndex = index;
									}
								index++;
							}
							results.println(unknown+" https://www.ncbi.nlm.nih.gov/pccompound/?term="+NISTMetabolite.getElementsByTagName("InChIKey").item(0).getTextContent().split(", ")[closestIndex]);
							
						}
						found = true;
					}
					
				}
				
				for(int temp = 0; temp < PUBCHEMList.size();temp++)
				{
					Element PUBCHEMMetabolite = PUBCHEMList.get(temp);
					Double PUBCHEMmz = Double.parseDouble(PUBCHEMMetabolite.getElementsByTagName("CenterMass").item(0).getTextContent().trim());
					if ((PUBCHEMmz-unknownmz)/unknownmz*1000000 < -PPMTolerance)
					{
						//do nothing
					}
					else if ((PUBCHEMmz-unknownmz)/unknownmz*1000000 > PPMTolerance)
					{
						//do nothing
					}
					else
					{
						String exactMass = PUBCHEMMetabolite.getElementsByTagName("ExactMass").item(0).getTextContent().trim();
						if(exactMass.indexOf(',')<0)
						{
							results.println(unknown+" https://www.ncbi.nlm.nih.gov/pccompound/?term="+PUBCHEMMetabolite.getElementsByTagName("InChIKey").item(0).getTextContent());
						}
						else
						{
							String[] parts = PUBCHEMMetabolite.getElementsByTagName("ExactMass").item(0).getTextContent().substring(1).split(", ");
							int index = 0;
							int closestIndex = 0;
							Double closestVal = 1.0;
							for(String clustermz : parts)
							{
								Double mzValue = Double.parseDouble(clustermz);
								if (mzValue - unknownmz < closestVal)
									{
										closestVal = mzValue;
										closestIndex = index;
									}
								index++;
							}
							results.println(unknown+" https://www.ncbi.nlm.nih.gov/pccompound/?term="+PUBCHEMMetabolite.getElementsByTagName("InChIKey").item(0).getTextContent().split(", ")[closestIndex]);
							
						}
						found = true;
					}
					
				}
				
				
				
				
				if(!found) //not matched to any database
				{
					//update results for NA
					results.println(unknown + ": was not matched to any of the databases.");
				}
			}
			
			results.close();
		} catch (Exception e) 
		{
			e.printStackTrace();
		}

	}

}
