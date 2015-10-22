package example;

public class Metabolite //created class to hold information relevant to PubChem metabolites
{
	private String Name;
	private String ID;
	private String Formula;
	private String Mass;
	private String InChIKey;
	public Metabolite(String Name, String ID, String Formula, String Mass, String InChIKey)
	{
		this.Name = Name;
		this.ID = ID;
		this.Formula = Formula;
		this.Mass = Mass;
		this.InChIKey = InChIKey;
	}
	public String getName()
	{
		return Name;
	}
	public String getID()
	{
		return ID;
	}
	public String getFormula()
	{
		return Formula;
	}
	public String getMass()
	{
		return Mass;
	}
	public String getInChIKey()
	{
		return InChIKey;
	}
	public String getCenterMass()
	{
		if(Mass.indexOf(',')<0)
		{
			return Mass;
		}
		else
		{
			String[] parts = Mass.split(",");
			Double firstVal = Double.parseDouble(parts[0]);
			Double lastVal = Double.parseDouble(parts[parts.length-1]);
			return Double.toString((firstVal+lastVal)/2);
		}
			
	}
}
