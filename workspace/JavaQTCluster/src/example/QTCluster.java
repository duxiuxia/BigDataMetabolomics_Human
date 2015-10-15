package example;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

import example.Metabolite;

public class QTCluster 
{
	public static void main(String[] args)
	{
		int PPMTolerance = 3;
		
		ArrayList<Metabolite> metaboliteList = new ArrayList<>();
		
		
		
		File dir = new File("/Users/MatthewDeitz/Desktop/MetaboliteDB");
		FilenameFilter filter = new FilenameFilter()
				{
					public boolean accept (File dir, String name) 
					{
						return name.endsWith(".csv");
					}
				};
		String[] str = dir.list(filter);
		for (String fileName : str)
		{
			BufferedReader fileReader = null;
			try
			{
				fileReader = new BufferedReader(new FileReader(fileName));
				String line = "";
				while ((line=fileReader.readLine())!=null)
				{
					String[] tokens = line.split(",");
					//check for NA or Null
					metaboliteList.add(new Metabolite(tokens[1],tokens[2],tokens[3],tokens[4],tokens[5]));
				}
				
			}
			catch (Exception e)
			{
				e.printStackTrace();
			}
			finally
			{
				try
				{
					fileReader.close();
				}
				catch(IOException e)
				{
					e.printStackTrace();
				}
			}
			
		}
		
		Double[] mzvalues = new Double[metaboliteList.size()];
		for(int temp = 0; temp < mzvalues.length;temp++)
		{
			Metabolite m = metaboliteList.get(temp);
			mzvalues[temp] = Double.parseDouble(m.getMass());
		}
		Arrays.sort(mzvalues);
		// need to clustermetaboliteList by mass
		
		/*
		 * 
		 * 
		 * methodtocontinue making clusters
		 * 
		 * 
		 */
		boolean flag = true;
		while(flag)
		{
			int index = findSmallestDist(mzvalues);
			if ((mzvalues[index+1]-mzvalues[index])/mzvalues[index]*1000000 > PPMTolerance)
			{
				flag = false;
			}
			else
			{
				
				int[] indicies = extendClusters(mzvalues,index,index+1, PPMTolerance);
				//print cluster to xml file
				
				//remove values from mz values and reset mz values
				ArrayList<Metabolite> cluster1 = new ArrayList<>();
				for(int x=indicies[0];x<=indicies[1];x++)
				{
					cluster1.add(metaboliteList.remove(x));
				}	
			}
		}
		
		
	}
	public static int[] extendClusters(Double[] mzvalues, int lowIndex, int highIndex, int PPMTolerance)
	{
		int[] indecies = {lowIndex,highIndex};
		
		
		if (lowIndex == 0)
		{
			Double centerMZ = (mzvalues[lowIndex] + mzvalues[highIndex+1])/2;
			if (centerMZ-mzvalues[lowIndex]/mzvalues[lowIndex]*1000000 > PPMTolerance)
			{
				return indecies;
			}
			else 
			{
				return extendClusters(mzvalues, lowIndex, highIndex+1,PPMTolerance);
			}
		}
		else if (highIndex +1 == mzvalues.length)
		{
			Double centerMZ = (mzvalues[lowIndex-1] + mzvalues[highIndex])/2;
			if (centerMZ-mzvalues[lowIndex-1]/mzvalues[lowIndex-1]*1000000 > PPMTolerance)
			{
				return indecies;
			}
			else 
			{
				return extendClusters(mzvalues, lowIndex-1, highIndex,PPMTolerance);
			}
		}
		else
		{
			Double centerMZ1 = (mzvalues[lowIndex-1] + mzvalues[highIndex])/2;
			Double centerMZ2 = (mzvalues[lowIndex] + mzvalues[highIndex+1])/2;
			Double centerMZ;
			int flag = 0;
			if (centerMZ1 < centerMZ2)
			{
				centerMZ = centerMZ1;
				flag = 1;
			}
			else
			{
				centerMZ = centerMZ2;
				flag = 2;
			}
			if (centerMZ-mzvalues[lowIndex-1]/mzvalues[lowIndex-1]*1000000 > PPMTolerance)
			{
				return indecies;
			}
			else 
			{
				if (flag ==1)
				{
					return extendClusters(mzvalues, lowIndex-1, highIndex,PPMTolerance);
				}
				else
				{
					return extendClusters(mzvalues, lowIndex, highIndex+1,PPMTolerance);
				}
				
			}
		}
	}
	public static int findSmallestDist(Double[] d)
	{
		Double smallestVal = d[1]-d[0];
		int index = 0;
		for (int temp = 0;temp < d.length-1;temp++)
		{
			if ((d[temp+1]-d[temp]) < smallestVal)
			{
				smallestVal = d[temp+1]-d[temp];
				index=temp;
			}
		}
		return index;
	}
}
