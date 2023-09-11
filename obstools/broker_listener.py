# -*- coding: utf-8 -*-
#
#  This file is part of LDTObserverTools.
#
#   This Source Code Form is subject to the terms of the Mozilla Public
#   License, v. 2.0. If a copy of the MPL was not distributed with this
#   file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
#  Created on 31-May-2023
#
#  @author: tbowers

"""LDTObserverTools contains python ports of various LDT Observer Tools

Lowell Discovery Telescope (Lowell Observatory: Flagstaff, AZ)
http://www.lowell.edu

This file contains the ``ActiveMQ_Listener`` class used by the
``deveny_collfocus`` routine for retrieving current values from the LDT's
ActiveMQ broker.  This is kept separate so that an ``ImportError`` on ``stomp``
simply causes this module to not import and the calling module can go about its
business.
"""

# Built-In Libraries

# 3rd-Party Libraries
import stomp
import yaml
import xmltodict

# Local Libraries


class ActiveMQ_Listener:
    """Broadcast (status) message listening class

    This class listens to the ActiveMQ broker and places a parsed version of
    the desired broker messages into attributes accessible by the calling
    function.  This class is adapted from similar ones used with the LORAX
    project, being developed at Lowell Observatory.

    Parameters
    ----------
    config_file : :obj:`str` or :obj:`pathlib.Path`
        The location of the configuration file to read
    """

    def __init__(self, config_file):
        # Initialize
        self.mounttemp_from_broker = {}
        self.grangle_from_broker = {}

        # Read the config file.
        with open(config_file, "r", encoding="utf-8") as stream:
            try:
                self.config = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)

        # Get the broker host from the configuration.
        # Make a connection to the broker.
        self.hosts = [tuple(self.config["broker_hosts"])]

        try:
            # Get a connection handle.
            self.conn = stomp.Connection(host_and_ports=self.hosts)

            # Set up a listener and and connect.
            self.conn.set_listener("", self.MyListener(self))
            self.conn.connect(wait=True)
        except stomp.exception.StompException:
            # If we cannot connect to the broker, return now and don't subscribe
            return

        # Subscribe to specified topics
        for topic in [
            f"{data_source}_incoming_topic" for data_source in ["mounttemp", "grangle"]
        ]:
            self.conn.subscribe(
                id=1,
                destination="/topic/" + self.config[topic],
                headers={},
            )

    class MyListener(stomp.ConnectionListener):
        """Class for listening to the STOMPing

        This listener runs asynchronously from the main body of the program
        through stomp's threading functionality.  When a broker message is
        received, it is parsed and placed into the ``parent`` attribute(s)
        to be accessed whenever the main code needs it.

        Parameters
        ----------
        parent : :obj:`class`
            The parent class into whose attributes values from the broker
            message are placed.
        """

        def __init__(self, parent):
            # Set parent for communication upward
            self.parent = parent

        def on_error(self, message):
            # TODO: This really needs to do something substantial
            print(f'received an error "{message}"')

        def on_message(self, message):
            # When a message is recieved, parse out and place into parent attributes

            # Get the topic
            topic = message.headers["destination"]

            if self.parent.config["mounttemp_incoming_topic"] in topic:
                self.parent.mounttemp_from_broker = {"MountTemp": float(message.body)}

            if self.parent.config["grangle_incoming_topic"] in topic:
                status = xmltodict.parse(message.body)["DevenyTelemetry"]
                self.parent.grangle_from_broker = self.parent.parse_deveny(status)

    def parse_deveny(self, status):
        """Parse an XMLTODICT status message

        Parameters
        ----------
        status : :obj:`dict`
            Translated XML status message

        Returns
        -------
        :obj:`dict`
            The parsed status dictionary for this particular use case.
        """
        # Return empty dictionary if no status
        if status is None:
            return {}

        ret_dict = dict(status.items())

        # Convert to float, if possible
        for key, val in status.items():
            try:
                ret_dict[key] = float(val)
            except ValueError:
                ret_dict[key] = val

        # Return
        return ret_dict
