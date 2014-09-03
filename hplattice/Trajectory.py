def open_file_stream(filename):
    """
    Helper function to open a file for writing.

    :param str filename: open this path for writing
    """
    return open(filename, 'w')


class Trajectory(object):
    """
    *Trajectory* objects write chain coordinates to output streams. They used
    to record the progress of an enumeration or monte carlo simulation.

    :param bool save_trajectory: ``True`` if trajectory should be saved to
                                 output stream.
    :param str trajectory_filename: write trajectory to this path
    :param open_stream_fcn: optional, call this function to open an output stream
    :type open_stream_fcn: callable
    """
    def __init__(self, save_trajectory, trajectory_filename,
                 open_stream_fcn=open_file_stream):
        self.save_trajectory = save_trajectory
        self.trajectory_filename = trajectory_filename
        if self.save_trajectory:
            self.output_stream = open_stream_fcn(self.trajectory_filename)
        else:
            self.output_stream = None
        self.frame_num = 0

    def snapshot(self, chain):
        """
        Record coordinates of chain to output stream.

        :param chain: Save coords of this chain.
        :type chain: :class:`hplattice.Chain.Chain`
        """
        if self.save_trajectory:
            # save chain coords in xyz format
            hp_string = chain.get_hp_string()
            coord_array = chain.get_coord_array()
            coord_strings = \
                ["%s\t%.1f\t%.1f\t%.1f\n" % (hp, row[0], row[1], 0.0) for hp, row in zip(hp_string, coord_array)]
            coord_str = "".join(coord_strings)
            num_atoms = coord_array.shape[0]
            self.output_stream.write("%d\nFrame %d\n" % (num_atoms, self.frame_num))
            self.output_stream.write(coord_str)
            self.frame_num += 1
        else:
            pass

    def finalize(self):
        """
        Close any open output streams.
        """
        if self.output_stream:
            self.output_stream.close()
        else:
            pass
